// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/hbonds/HBondSet.hh>

#include <core/sequence/util.hh>

//Mmmm.. constraints.
//#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/PeriodicFunc.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
//#include <core/scoring/constraints/AngleConstraint.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/NonIdealProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>
#include <protocols/idealize.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


#include <core/scoring/TenANeighborGraph.hh>

//Backrub
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <core/scoring/dssp/Dssp.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/OptionKeys.OptionKeys.gen.hh>
#include <core/options/keys/run.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/backrub.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Following was developed in 2008 for CASP8 by R. Das for targets with unambiguous, nearly continuous
// templates.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, chi_stats )
OPT_KEY( Boolean, copy_native_chi )
OPT_KEY( Boolean, rhiju_fold_tree )
OPT_KEY( Boolean, use_native_CA )
OPT_KEY( Boolean, vary_geometry )
OPT_KEY( Boolean, mask_loop )
OPT_KEY( Boolean, backrub_test )
OPT_KEY( Boolean, easy_loop_model )
OPT_KEY( Boolean, cst_relax )
OPT_KEY( Boolean, t469_zinc_tether )
OPT_KEY( String, soft_segment_file )
OPT_KEY( String, sequence_mask_file )
OPT_KEY( String, secstruct_file )
OPT_KEY( String, cst_file )
OPT_KEY( Real, chi1_constraint_weight )
OPT_KEY( Real, chi2_constraint_weight )
OPT_KEY( Real , CA_tether )
OPT_KEY( Real , soft_CA_tether )
OPT_KEY( Integer, rounds )
OPT_KEY( Integer, cst_trim_loop )

OPT_1GRP_KEY(IntegerVector, backrub, pivot_residues)
OPT_1GRP_KEY(StringVector, backrub, pivot_atoms)
OPT_1GRP_KEY(Integer, backrub, min_atoms)
OPT_1GRP_KEY(Integer, backrub, max_atoms)
OPT_1GRP_KEY(Integer, backrub, ntrials)
OPT_1GRP_KEY(Real, backrub, sc_prob)
OPT_1GRP_KEY(Real, backrub, sc_prob_uniform)
OPT_1GRP_KEY(Real, backrub, mc_kt)
OPT_1GRP_KEY(Real, backrub, mm_bend_weight)


static numeric::random::RandomGenerator RG(22257);

//////////////////////////////////////////////////////////////////
Size
read_alignment_fasta_file(
				utility::vector1< std::string >  & sequences,
				utility::vector1< std::string >  & pdb_names,
				std::string const & fasta_file
 )
{
	sequences.clear();
	pdb_names.clear();

	utility::io::izstream data_stream( fasta_file );

	if ( data_stream.fail() )		utility_exit_with_message( "Could not find fasta ifile: " + fasta_file  );


	std::string line;
	std::string sequence = "";

	while( getline(data_stream, line) 	) {

		while ( line.size() == 0 || line[0] == '>' || line[0] == ' '  ) {

			//Finished a sequence?
			if ( sequence.size() > 0 ) 	sequences.push_back( sequence );

			//Is this a tag with the pdb names?
			if ( line.size() > 0 && line[0] == '>' ) {
				pdb_names.push_back( line.substr( 1, line.size()-1 ) );
			}

			sequence = "";
			getline(data_stream, line);
			if ( data_stream.fail() ) break;
		}

		if ( data_stream.fail() ) break;

		//We're inside a sequence.
		sequence += line;
	}
	if (sequence.size() > 0 ) sequences.push_back( sequence );

	if ( sequences.size() == 0 )		utility_exit_with_message( "Could not find sequence inside file" + fasta_file  );

	Size const alignment_length = sequences[1].size();
	for (Size n = 1; n <= sequences.size(); n++ ){
		std::cout << pdb_names[n] << "  " << sequences[n] << std::endl;
		assert( sequences[n].size() == 	alignment_length );
	}

	return alignment_length;
}

///////////////////////////////////////////////////////////////////////
void
setup_mask(
  FArray1D_bool & sequence_mask,
	utility::vector1< std::string > const & sequences )
{

	sequence_mask = false;

	Size const alignment_length = sequences[1].size();

	///////////////////////////////////////////////////////////////
	// First pass -- rule out any part that is gapped in the alignment.
	for (Size i = 1; i <= alignment_length; i++ ) {

		bool found_a_gap( false );
		for (Size n = 1; n <= sequences.size(); n++ ){
			if ( sequences[n][i-1] == '-' ) {
				found_a_gap = true;
				break;
			}
		}

		if (!found_a_gap) {
			sequence_mask( i ) = true;
		}

	}

	///////////////////////////////////////////////////////////////
	// Second pass -- look for ungapped part that are bracketed by gapped regions
	Size const look_for_gap( 4 );
	for (int i = 1; i <= int(alignment_length); i++ ) {

		if ( sequence_mask(i) == false ) continue; //Don't worry about it.

		bool found_gap_before( false ), found_gap_after( false );
		for (int offset = -1 * look_for_gap; offset < 0; offset ++ ) {
			if ( i+offset > 1 && sequence_mask( i+offset ) == false ) {
				found_gap_before = true;
				break;
			}
		}
		if (!found_gap_before) continue;

		for (int offset = 1; Size(offset) <= look_for_gap; offset ++ ) {
			if ( i+offset <= int(alignment_length) && sequence_mask( i+offset ) == false ) {
				found_gap_after = true;
				break;
			}
		}
		if (!found_gap_after) continue;

		std::cout << "MASK: Region that is not nominally gapped but looks fishy: " << i << std::endl;
		sequence_mask( i ) = false;
	}

	////////////////////////////////////////////////
	// Debug output
	std::cout << "SETUP MASK ==>  exclude ";
	for (Size i = 1; i <= alignment_length; i++ ) {
		if ( sequence_mask(i) == false ) std::cout << i << " " ;
	}
	std::cout << std::endl;

	////////////////////////////////////////////////
	// Save into a file.
	Size count( 0 );
	utility::io::ozstream out( options::option[ options::OptionKeys::sequence_mask_file ]  );
	for (Size i = 1; i <= alignment_length; i++ ) {
		if (sequences[1][i-1] == '-' ) continue; //Assume native numbering
		count++;
		out << sequence_mask(i);
	}
	out << std::endl;
	out.close();

}

////////////////////////////////////////////////////////////////////////////
void
setup_alignment_map( utility::vector1< std::map< Size, Size > > & alignment2sequence, utility::vector1 < std::string > const & sequences )
{

	//Figure out mapping of alignment to each sequence.

	Size const alignment_length( sequences[1].size() );
	for (Size n = 1; n <= sequences.size(); n++ ){
		Size count( 0 );
		std::map< Size, Size > mapping;
		for (Size i = 1; i <= alignment_length; i++ ) {
			if ( sequences[n][i-1] != '-' ) {
				count++;
				mapping[ i ] = count;
			} else {
				mapping[ i ] = 0;
			}
		}
		//		std::cout << "TOTAL RES " << n << " --> " << count << std::endl;
		alignment2sequence.push_back( mapping );
	}

}


////////////////////////////////////////////////////////////////////////////////////
std::string
apply_mask(
					 FArray1D_bool & sequence_mask,
					 utility::vector1< std::map< Size, Size > > & alignment2sequence,
					 utility::vector1< std::string > const & pdb_names,
					 std::string const & which_file,
					 pose::Pose & pose )
{

	pose::Pose temp_pose;

	temp_pose.clear();

	// Find which alignment to use, based on input pdb name.
	Size which_sequence( 0 );
	for( Size n=1; n<=pdb_names.size(); n++ ){
		//		std::cout << pdb_names[n] << " " << which_file << std::endl;
		if ( pdb_names[n] == which_file ) {
			which_sequence = n;
			break;
		}
	}

	if ( which_sequence == 0 ) utility_exit_with_message( "Problem with finding tag " + which_file + " in fasta file." );


	// Create starting pdb, apply mask.
	Size const alignment_length( alignment2sequence[1].size() );
	bool in_cutpoint( false );
	Size count( 0 );
	for (Size i = 1; i <= alignment_length; i++ ){
		if (  sequence_mask( i )  ) {
			Size const pdb_number( alignment2sequence[ which_sequence ][i] );
			count++;
			//			std::cout << "MASK " << i << " " << sequence_mask(i) << " " << pdb_number <<  std::endl;
			if (in_cutpoint) {
				temp_pose.append_residue_by_jump( pose.residue( pdb_number ), count-1 );
			} else {
				temp_pose.append_residue_by_bond( pose.residue( pdb_number ) );
			}

			temp_pose.set_secstruct( count, pose.secstruct( pdb_number ) );
			in_cutpoint = false;

		} else {
			if (count > 0 ) in_cutpoint = true;
		}
	}

	std::cout << "FOLD TREE" << temp_pose.fold_tree() << std::endl;
	pose = temp_pose;

	std::string masked_sequence = pose.sequence();
	return masked_sequence;

}

/////////////////////////////////////////////////////////////////////////////////
// This could easily be alignized...
void
setup_t469_zn_tether( pose::Pose & pose ){
	using namespace id;
	using namespace scoring::constraints;
	using namespace options;
	using namespace options::OptionKeys;

	ConstraintSetOP cst_set( 	pose.constraint_set()->clone() ) ;

	Real const tetrahedron_optimum_distance( 4.0);
	Real const tetrahedron_sd( 1.0 ); //reasonably loose.
	FuncOP tether_func( new HarmonicFunc( tetrahedron_optimum_distance, tetrahedron_sd ) );

	Size const his1( 48 ), cys2( 30 ), cys3( 33 ), glu4( 39 );

	// Don't actually know which nitrogen on His. Hmm
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(his1).atom_index( " NE2" ), his1 ),
															AtomID( pose.residue(cys2).atom_index( " SG " ), cys2 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(his1).atom_index( " NE2" ), his1 ),
															AtomID( pose.residue(cys3).atom_index( " SG " ), cys3 ),
															tether_func ) );
	// Don't know why oxygen on glu.
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(his1).atom_index( " NE2" ), his1 ),
															AtomID( pose.residue(glu4).atom_index( " OE1" ), glu4 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(his1).atom_index( " NE2" ), his1 ),
															AtomID( pose.residue(glu4).atom_index( " OE2" ), glu4 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(cys2).atom_index( " SG " ), cys2 ),
															AtomID( pose.residue(cys3).atom_index( " SG " ), cys3 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(cys2).atom_index( " SG " ), cys2 ),
															AtomID( pose.residue(glu4).atom_index( " OE1" ), glu4 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(cys2).atom_index( " SG " ), cys2 ),
															AtomID( pose.residue(glu4).atom_index( " OE2" ), glu4 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(cys3).atom_index( " SG " ), cys3 ),
															AtomID( pose.residue(glu4).atom_index( " OE1" ), glu4 ),
															tether_func ) );
	cst_set->add_constraint(
			new AtomPairConstraint( AtomID( pose.residue(cys3).atom_index( " SG " ), cys3 ),
															AtomID( pose.residue(glu4).atom_index( " OE2" ), glu4 ),
															tether_func ) );

	pose.constraint_set( cst_set );
}


////////////////////////////////////////////////////////////////////////////
void
prepare_start_model( pose::Pose const & template_pose,
										 pose::Pose const & pose_with_desired_sequence,
										 scoring::ScoreFunction & scorefxn_input,
										 pose::Pose & pose,
										 FArray1D_bool & conserved )
{

	using namespace scoring;
	using namespace scoring::constraints;
	using namespace options;
	using namespace options::OptionKeys;

	scoring::ScoreFunction scorefxn = scorefxn_input;

	pose = template_pose;
	conserved = false;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {

		utility::vector1< bool > repack_aalist( chemical::num_canonical_aas, false );
		repack_aalist[  pose_with_desired_sequence.residue(i).aa() ] = true;
	 	task->nonconst_residue_task(i).restrict_absent_canonical_aas( repack_aalist );

		//Disallow any change in conserved residues.
		if ( pose.residue(i).aa() == pose_with_desired_sequence.residue(i).aa() ) {
			//			std::cout << "CONSERVED RESIDUE: " << i << std::endl;
			conserved( i ) = true;

			// Hold on -- this could be a big problem since atom ids change if residue changes.
			//Put on a constraint to lower the probability of this chi1 changing.
	// 		if ( pose.residue(i).nchi() > 0 && pose_with_desired_sequence.residue(i).nchi()) {
// 				id::TorsionID my_ID( i, id::CHI, 1 );
// 				id::AtomID id1,id2,id3,id4;
// 				pose.conformation().get_torsion_angle_atom_ids( my_ID, id1, id2, id3, id4 );
// 				pose.constraint_set()->add_constraint(
// 																							new DihedralConstraint( id1, id2, id3, id4,
// 																																			new PeriodicFunc( numeric::conversions::radians( pose.torsion( my_ID ) ),
// 																																												option[ chi1_constraint_weight ], 1 ) ) );
// 				scorefxn.set_weight( dihedral_constraint, 1.0  );
// 			}

			//task->nonconst_residue_task(i).prevent_repacking();
		}



//  			//Put on a constraint to lower the probability of this chi2 changing.
// 			if ( pose.residue(i).nchi() > 1 ) {
// 				id::TorsionID my_ID( i, id::CHI, 2 );
// 				id::AtomID id1,id2,id3,id4;
// 				pose.conformation().get_torsion_angle_atom_ids( my_ID, id1, id2, id3, id4 );
// 				pose.constraint_set()->add_constraint(
//   							 new DihedralConstraint( id1, id2, id3, id4,
//   							 new PeriodicFunc( numeric::conversions::radians( pose.torsion( my_ID ) ),
//   																 option[ chi2_constraint_weight ], 1 ) ) );
// 			}

	}


	//	scorefxn.set_weight( fa_rep, 0.01 );

	if ( option[ t469_zinc_tether]() ) setup_t469_zn_tether( pose );

	pack::pack_rotamers( pose, scorefxn, task);

	//Copy chi's from native... best case scenario
	if ( option[ copy_native_chi ]() ){
		for (Size i = 1; i <= pose.total_residue(); i++) {
			for ( Size n = 1; n <= pose.residue(i).nchi(); n++ ){
				pose.set_chi( n, i, pose_with_desired_sequence.chi( n, i ) ) ;
			}
		}
	} else {
	}

}



////////////////////////////////////////////////////////////////////////////
void
prepare_start_model( pose::Pose const & template_pose,
										 std::string desired_sequence,
										 scoring::ScoreFunction & scorefxn_input,
										 pose::Pose & pose,
										 FArray1D_bool & conserved )
{

	using namespace scoring;
	using namespace scoring::constraints;
	using namespace options;
	using namespace options::OptionKeys;

	scoring::ScoreFunction scorefxn = scorefxn_input;

	pose = template_pose;
	conserved = false;

	std::cout << "HEY!! " << template_pose.sequence() << std::endl;
	std::cout << "HEY!! " << desired_sequence << std::endl;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		utility::vector1< bool > repack_aalist( chemical::num_canonical_aas, false );
		repack_aalist[  core::chemical::aa_from_oneletter_code( desired_sequence[i-1] ) ] = true;
	 	task->nonconst_residue_task(i).restrict_absent_canonical_aas( repack_aalist );
	}

	pack::pack_rotamers( pose, scorefxn, task);


}

////////////////////////////////////////////////////////////////////////////
void
repack(	 pose::Pose & pose, scoring::ScoreFunction const & scorefxn )
{
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace options;
	using namespace options::OptionKeys;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->disallow_quench( true );
	for (Size i = 1; i <= pose.total_residue(); i++) {
	 	task->nonconst_residue_task(i).restrict_to_repacking();
	}

	pack::pack_rotamers( pose, scorefxn, task);

}

////////////////////////////////////////////////////////////////////////////
void
check_chi_correct( pose::Pose const & native_pose,
									 pose::Pose const & pose,
									 Size const & i,
									 bool & chi1_exists,
									 bool & chi2_exists,
									 bool & chi1_correct,
									 bool & chi1_chi2_correct )
{

	static Real const CHI_CUTOFF( 30.0 );

	chi1_exists = false;
	chi2_exists = false;
	chi1_correct = false;
	chi1_chi2_correct = false;

	if ( pose.residue(i).nchi() > 0 ) {
		chi1_exists = true;
		if ( std::abs( native_pose.chi( 1, i ) - pose.chi( 1, i ) ) < CHI_CUTOFF ) chi1_correct = true;

		if ( pose.residue(i).nchi() > 1 ) {
			chi2_exists = true;
			if ( chi1_correct && std::abs( native_pose.chi( 2, i ) - pose.chi( 2, i ) ) < CHI_CUTOFF ) chi1_chi2_correct = true;
		}

	}



}

////////////////////////////////////////////////////////////////////////////
void
output_backbone_stats( pose::Pose & native_pose, pose::Pose & pose )
{

	utility::io::ozstream out( "backbone_stats.txt"  );
	Size const nres( native_pose.total_residue() );
	for (Size i = 1; i < nres; i++ ){
		out << "BB_STATS " <<
			" " << native_pose.phi(i) << " " << pose.phi(i) <<
			" " << native_pose.psi(i) << " " << pose.psi(i) <<
			" " << native_pose.omega(i) << " " << pose.omega(i) << std::endl;
	}

}

////////////////////////////////////////////////////////////////////////////
void
output_chi_stats( pose::Pose & native_pose, pose::Pose & pose,
									FArray1D_bool sequence_mask,
									utility::io::ozstream & out,
									std::string const & pdb_file, Size const n)
{
	Size const pair_score_cb_thresh( 16 );

	Size const nres( native_pose.total_residue() );

	Size n_chi1(0), n_chi1_buried(0), n_chi1_exposed(0);
	Size n_chi1_correct(0), n_chi1_correct_buried(0), n_chi1_correct_exposed(0);
	Size n_chi2(0), n_chi2_buried(0), n_chi2_exposed(0);
	Size n_chi1_chi2_correct(0), n_chi1_chi2_correct_buried(0), n_chi1_chi2_correct_exposed(0);

	if (n == 1 ) {
		out << A( 30, "PDB" ) << " ==> "
				<< A( 8, "chi1_bur" ) << " "
				<< A( 8, "chi1_exp" ) << " "
				<< A( 8, "chi1" ) << "    "
				<< A( 8, "chi2_bur" ) << " "
				<< A( 8, "chi2_exp" ) << " "
				<< A( 8, "chi2" ) << std::endl;
	}

	//Check burial.
	scoring::TenANeighborGraph const & neighbor_graph
		( native_pose.energies().tenA_neighbor_graph() );

	//This doesn't do LAST RESIDUE!! Cause of issues in last residues in server models?
	bool problem( false );

	//For Rasmol visualization...
	bool found_wrong = false;
	std::cout << "CHI1 WRONG FOR BURIED RESIDUE: ";

	for (Size i = 1; i < nres; i++ ){

		if  ( !sequence_mask(i) ) continue;

		if ( pose.sequence()[i-1] != native_pose.sequence()[i-1] ) {
			problem = true;
			break;
		}

		Size const nb = neighbor_graph.get_node( i )->num_neighbors_counting_self();

		bool buried( false );
		if ( nb > pair_score_cb_thresh ) buried = true;

		bool chi1_exists( false ), chi2_exists( false ), chi1_correct( false ), chi1_chi2_correct( false );

		check_chi_correct( native_pose, pose, i,
											 chi1_exists, chi2_exists,
											 chi1_correct, chi1_chi2_correct );

		if ( chi1_exists ) {
			n_chi1++;
			if (chi1_correct) n_chi1_correct++;

			if (buried) {
				n_chi1_buried++;
				if (chi1_correct) n_chi1_correct_buried++;

				if (!chi1_correct) {
					if (found_wrong) std::cout << ",";
					found_wrong = true;
					std::cout << i;
				}

			} else {
				n_chi1_exposed++;
				if (chi1_correct) n_chi1_correct_exposed++;
			}
		}

		if ( chi2_exists ) {
			n_chi2++;
			if (chi1_chi2_correct) n_chi1_chi2_correct++;

			if (buried) {
				n_chi2_buried++;
				if (chi1_chi2_correct) n_chi1_chi2_correct_buried++;
			} else {
				n_chi2_exposed++;
				if (chi1_chi2_correct) n_chi1_chi2_correct_exposed++;
			}
		}

	}

	std::cout << std::endl;

	if (problem) return;

	out << A( 30, pdb_file ) << " ==> " ;
	out
		<< F(8,4, n_chi1_correct_buried/(1.0 * n_chi1_buried) ) << " "
		<< F(8,4, n_chi1_correct_exposed/(1.0 * n_chi1_exposed) ) << " "
		<< F(8,4, n_chi1_correct/(1.0 * n_chi1) ) << "    "
		<< F(8,4, n_chi1_chi2_correct_buried/(1.0 * n_chi2_buried) ) << " "
		<< F(8,4, n_chi1_chi2_correct_exposed/(1.0 * n_chi2_exposed) ) << " "
		<< F(8,4, n_chi1_chi2_correct/(1.0 * n_chi2) ) << std::endl;

}


////////////////////////////////////////////////////////////////////////////
void
setup_secstruct( pose::Pose & template_pose, std::string const & template_secstruct_file )
{
	utility::io::izstream data_stream( template_secstruct_file );
	std::string line;
	getline(data_stream, line);
	for (Size i = 1; i <= template_pose.total_residue(); i++ ) {
		template_pose.set_secstruct(i,  line[i-1] );
	}
}

////////////////////////////////////////////////////////////////////////////
void
setup_secstruct_dssp( pose::Pose & pose )
{
	core::scoring::dssp::Dssp dssp( pose );
	FArray1D_char dssp_secstruct( pose.total_residue() );
	dssp.dssp_reduced( dssp_secstruct );
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		pose.set_secstruct(i,  dssp_secstruct(i) );
	}

}


////////////////////////////////////////////////////////////////////////////
Size
find_middle_of_vector( utility::vector1 <Size> v){
	Size value = v[ Size( 0.5 * v.size() + 0.5 ) ];
	//	std::cout << "MIDDLE " << value << std::endl;
	return value;
}

////////////////////////////////////////////////////////////////////////////
void
setup_rhiju_fold_tree( pose::Pose & pose )
{

	Size const nres( pose.total_residue() );

	//Quick check
	std::cout <<  "SECSTRUCT: " ;
	for (Size i = 1; i <= nres; i++ ) {
		std::cout << pose.secstruct( i );
	}
	std::cout << std::endl;

	// Make list of good solid stretches of secondary structures.
	FArray1D_bool good_secstruct( nres,  false );
	for (Size i = 1; i <= nres; i++ ) {
		if ( pose.secstruct(i) != 'L' ) good_secstruct( i ) = true;
	}

	// Remove loners.
	for (Size i = 2; i < nres; i++ ) {
		if ( pose.secstruct(i-1) == 'L' && pose.secstruct(i+1)=='L' ) {
			good_secstruct( i ) = false;
		}
	}

	utility::vector1 < std::pair< Size, Size> > secstruct_chunks;
	Size start( 0 ), end( 0 );
	bool in_secstruct( false );
	for (Size i = 1; i <= nres; i++ ) {
		if ( good_secstruct(i) && !in_secstruct  ){
			in_secstruct = true;
			start = i;
		}
		if ( (!good_secstruct(i)  || pose.fold_tree().is_cutpoint( i )) && in_secstruct  ){
			in_secstruct = false;
			end = i;
			secstruct_chunks.push_back( std::make_pair( start, end ) );
		}
	}

	Size const num_segments( secstruct_chunks.size()  );
	std::cout << "Found " << num_segments << " secondary structure segments." << std::endl;

	// Make list of jumps (which will go to an appended virtual residue)
	// and potential cutpoints.
	Size num_cutpoints( num_segments - 1 );
	ObjexxFCL::FArray2D <int> jump_points( 2, num_cutpoints );
	ObjexxFCL::FArray1D <int> cuts( num_cutpoints );
	FArray1D_bool obligate_cutpoint( num_cutpoints, false );

	bool const star_fold_tree( false );

	if ( star_fold_tree ) {
		num_cutpoints = num_segments;
		jump_points.dimension( 2, num_cutpoints );
		cuts.dimension( num_cutpoints );
		obligate_cutpoint.dimension( num_cutpoints, false );

		// Make list of jumps (which will go to an appended virtual residue)
		// and potential cutpoints.
		ObjexxFCL::FArray2D <int> jump_points( 2, num_segments );
		ObjexxFCL::FArray1D <int> cuts( num_segments );
		FArray1D_bool obligate_cutpoint( num_segments, false );
		for (Size n = 1; n <= num_segments; n++ ) {
			jump_points( 1, n ) = ( secstruct_chunks[n].first + secstruct_chunks[n].second ) / 2;
			jump_points( 2, n ) = nres + 1;

			if ( n < num_segments ) {
				Size found_cutpoint( 0 );
				for (Size i = secstruct_chunks[ n ].second; i < secstruct_chunks[ n+1 ].first; i++ ) {
					if ( pose.fold_tree().is_cutpoint( i ) ) {
						found_cutpoint = i;
						obligate_cutpoint( n ) = true;
						break;
					}
				}
				if (found_cutpoint > 0)  {
					cuts( n ) = found_cutpoint;
			} else {
					cuts( n ) = ( secstruct_chunks[ n ].second + secstruct_chunks[ n+1 ].first ) / 2;
				}
			} else {
				cuts( n ) = nres;
			}
		}

	} else { //pretty straightforward fold tree with internal connection points.

		for (Size n = 1; n <= num_segments; n++ ) {

			Size found_cutpoint( 0 );
			for (Size i = secstruct_chunks[ n ].second; i < secstruct_chunks[ n+1 ].first; i++ ) {
				if ( pose.fold_tree().is_cutpoint( i ) ) {
					found_cutpoint = i;
					obligate_cutpoint( n ) = true;
					break;
				}
			}
			if (found_cutpoint > 0)  {
				cuts( n ) = found_cutpoint;
			} else {
				cuts( n ) = ( secstruct_chunks[ n ].second + secstruct_chunks[ n+1 ].first ) / 2;
			}
		}

		//Find potential jump points.
		//Start by defining residues making contacts.
		//  Pairs of segments with largest interfaces (most residues in contact)
		//   will be favored to be connected by jumps.
		//  Jumps will connect residues that are "in the middle" of each interface.
		bool const verbose = true;

		FArray2D_int num_contact_pairs( num_segments, num_segments, 0);
		FArray3D_int potential_jump_points( 2, num_segments, num_segments, 0);
		Real const CUTOFF_DISTANCE2( 8.0*8.0 );

		for (Size i = 1; i <= num_segments; i++ ){
			for (Size j = i+1; j <= num_segments; j++ ){

				Size temp_num_contact_pairs( 0 );
				utility::vector1< Size > interface_points_i;
				utility::vector1< Size > interface_points_j;

				//Just to avoid problems with extended structures, start out with 1 "contact pair"
				// between the middles of the two segments.
				temp_num_contact_pairs++;
				interface_points_i.push_back(  (secstruct_chunks[i].first + secstruct_chunks[i].second)/2 );
				interface_points_j.push_back(  (secstruct_chunks[j].first + secstruct_chunks[j].second)/2 );

				for (Size pos1 = secstruct_chunks[i].first; pos1 <= secstruct_chunks[i].second; pos1++ ) {
					for (Size pos2 = secstruct_chunks[j].first; pos2 <= secstruct_chunks[j].second; pos2++ ) {
						Real const distance2  = ( pose.residue(pos1).nbr_atom_xyz() - pose.residue(pos2).nbr_atom_xyz() ).length_squared();
						if (distance2  < CUTOFF_DISTANCE2 ){
							temp_num_contact_pairs++;
							interface_points_i.push_back( pos1 );
							interface_points_j.push_back( pos2 );
						}
					} // pos2
				} // pos1

				//How we doing? Keep things symmetrized.
				num_contact_pairs( i, j) = temp_num_contact_pairs;
				num_contact_pairs( j, i) = temp_num_contact_pairs;

				potential_jump_points(1, i, j) = find_middle_of_vector( interface_points_i );
				potential_jump_points(2, i, j) = find_middle_of_vector( interface_points_j );
				//upstream jump point > downstream jump point.
				potential_jump_points(1, j, i) = find_middle_of_vector( interface_points_i );
				potential_jump_points(2, j, i) = find_middle_of_vector( interface_points_j );

				if (verbose) std::cout << "Number of contact pairs " << i << " and " << j << " : " << num_contact_pairs(i,j) <<std::endl;
			} // j
		} // i

		//Sort-of neighbor joining...
		FArray1D_int segments_so_far( num_segments, 0);
		FArray1D_bool segment_in_tree( num_segments, false);

		//To start the tree, pick the two segments with the most extesnive interface.
		int num_contact_pairs_best = 0;
		int connection_point( 0 ), new_segment( 0 );
		for (Size i = 1; i <= num_segments; i++ ){
			for (Size j = i+1; j <= num_segments; j++ ){
				if (num_contact_pairs(i,j) > num_contact_pairs_best){
					num_contact_pairs_best = num_contact_pairs(i,j);
					connection_point = i;
					new_segment = j;
				}
			}
		}
		segments_so_far(1) = connection_point;
		segments_so_far(2) = new_segment;
		segment_in_tree( connection_point ) = true;
		segment_in_tree( new_segment ) = true;
		jump_points( 1, 1) = potential_jump_points( 1, connection_point, new_segment);
		jump_points( 2, 1) = potential_jump_points( 2, connection_point, new_segment);;

		//Now agglomerate the remaining segments into this tree, one-by-one
		for (Size count = 2; count<=num_cutpoints; count++){
			Size num_contact_pairs_best = 0;
			Size connection_point = 0;
			Size new_segment = 0;
			for (Size i = 1; i <= count; i++){
				int const current_segment = segments_so_far(i);
				for (Size j = 1; j <= num_segments; j++){
					if (segment_in_tree(j)) continue;
					if ( Size( num_contact_pairs(i,j) ) > num_contact_pairs_best){
						num_contact_pairs_best = num_contact_pairs(i,j);
						new_segment = j;
						connection_point = current_segment;
					}
				}
			}
			segments_so_far( count+1 ) = new_segment;
			segment_in_tree( new_segment ) = true;
			jump_points( 1, count ) = potential_jump_points( 1, connection_point, new_segment );
			jump_points( 2, count ) = potential_jump_points( 2, connection_point, new_segment );;
		}


		//Need to reorder jump_points to be in sequential order, or tree_from_jumps_and_cuts
		// gets confused.
		for (Size count = 1; count <= num_cutpoints; count++ ){
			if (jump_points(1, count) > jump_points(2, count) ) {
				int tmp = jump_points(1, count);
				jump_points(1, count) = jump_points(2, count);
				jump_points(2, count) = tmp;
			}
			if (verbose) std::cout << "Want a jump from " << jump_points(1, count) << " to " << jump_points(2, count) << std::endl;
		}

	}

	pose::Pose start_pose;
	start_pose = pose;

	// Create new fold tree. Try the virtual residue trick?
	using namespace conformation;
	if ( star_fold_tree ) {
		ResidueOP rsd( ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ) );
		pose.append_residue_by_jump( *rsd, 1 );
	}

	kinematics::FoldTree f( pose.total_residue() );
	f.tree_from_jumps_and_cuts( pose.total_residue(), num_cutpoints,
															jump_points, cuts, 1, false /*verbose*/ );
	pose.fold_tree( f );

	std::cout << "NEW! " << pose.fold_tree() << std::endl;

	//OK, let's do the chainbreak variant thing, except at obligate cutpoints.
	for (Size n = 1; n < num_segments; n++ ){
		if ( obligate_cutpoint( n ) ) continue;
		chemical::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cuts(n) );
		chemical::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cuts(n)+1 );
		pose.set_psi  ( cuts(n),   start_pose.psi( cuts(n) ) );
		pose.set_omega( cuts(n),   start_pose.omega( cuts(n) ) );
		pose.set_phi  ( cuts(n+1), start_pose.phi( cuts(n+1) ) );
	}



}

///////////////////////////////////////////////////////////////
void
read_soft_segment_file( FArray1D_bool & in_soft_segment ) {

	using namespace core::options;
	using namespace core::options::OptionKeys;

	in_soft_segment = false;
	if ( !option[ soft_segment_file ].user() ) return;

	utility::io::izstream data_stream( option[ soft_segment_file]() );
	std::string line;
	getline(data_stream, line);
	std::istringstream line_stream( line );
	Size pos1, pos2;
		while( !line_stream.fail() ) {
			line_stream >> pos1 >> pos2;
			for (Size i = pos1; i <= pos2; i++ ) in_soft_segment(i) = true;
		}

}

////////////////////////////////////////////////////////////////////////////
void
setup_CA_constraints( pose::Pose & pose, pose::Pose const & src_pose ) {

	using namespace id;
	using namespace scoring::constraints;
	using namespace options;
	using namespace options::OptionKeys;

	static Real const CA_cutoff( 9.0 );

	Size const nres( src_pose.total_residue() );

	//Don't constrain residues within 2 residues of a gap/cutpoint...
	// similar to loop-rebuild protocol
	FArray1D_bool allow_constraint( nres, true );
	static int const DISTANCE_TO_GAP_CUTOFF( options::option[ options::OptionKeys::cst_trim_loop ] );
	for (int i = 1; i <= int(nres); i++ ) {

		for (int offset = 0; offset < DISTANCE_TO_GAP_CUTOFF; offset++ )  {
			if ( ( (i+offset) <= int(nres)) && src_pose.fold_tree().is_cutpoint( i+offset) ) allow_constraint( i ) = false;
		}
		for (int offset = -DISTANCE_TO_GAP_CUTOFF; i < 0; offset++ )  {
			if ( ( (i+offset) >=1)  && src_pose.fold_tree().is_cutpoint( i+offset) ) allow_constraint( i ) = false;
		}
		if ( ! allow_constraint( i ) ) {
			std::cout << "NOT PUTTING CONSTRAINTS ON RESIDUE: " << i << std::endl;
		}
	}


	//User can also specify segments that shouldn't be tethered as strongly
	FArray1D_bool in_soft_segment( nres, false );
	read_soft_segment_file( in_soft_segment );

	static Real const HARD_CONSTRAINT_SD = option[ CA_tether ]();
	static Real const SOFT_CONSTRAINT_SD = option[ soft_CA_tether]();

	ConstraintSetOP cst_set( 	pose.constraint_set()->clone() ) ;

	for (Size i = 1; i <= nres; i++ ) {

		Vector const CA_i( src_pose.residue(i).xyz( " CA " ) );

		if ( !allow_constraint( i ) ) continue;

		for (Size j = i+1; j <= nres; j++ ) {

			if ( !allow_constraint( j ) ) continue;

			Vector const CA_j( src_pose.residue(j).xyz( " CA " ) );

			Real const CA_dist = (CA_i - CA_j).length();
			if ( CA_dist < CA_cutoff ) {

					if ( in_soft_segment(i) || in_soft_segment(j) )  {
						//						std::cout << "SOFT?" << i << " " << j << std::endl;
						//continue;
						cst_set->add_constraint(
																									new AtomPairConstraint( AtomID( pose.residue(i).atom_index( " CA " ), i ),
																																					AtomID( pose.residue(j).atom_index( " CA " ), j ),
																																					new HarmonicFunc( CA_dist, SOFT_CONSTRAINT_SD ) ) );
					} else {
						cst_set->add_constraint(
																									new AtomPairConstraint( AtomID( pose.residue(i).atom_index( " CA " ), i ),
																																					AtomID( pose.residue(j).atom_index( " CA " ), j ),
																																					new HarmonicFunc( CA_dist, HARD_CONSTRAINT_SD ) ) );
					}
			}

		}

	}

	pose.constraint_set( cst_set );

}


/////////////////////////////////////////////////////////////////////////////////
void vary_geometry_sidechains( pose::Pose & pose, kinematics::MoveMap & mm )
{

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set( 	pose.constraint_set()->clone() ) ;

	//Change this to also include D DOF's for sidechains.
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		for (Size n = 1; n <= pose.residue(i).nchi(); n++ ) {
			TorsionID tor_id( i, CHI, n );
			AtomID id1,id2,id3,id4, my_ID;
			pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );

			core::kinematics::tree::AtomCOP atom2 ( & pose.atom_tree().atom( id2 ) );
			core::kinematics::tree::AtomCOP atom3 ( & pose.atom_tree().atom( id3 ) );

			DOF_ID dof_id;
			if ( atom2->parent() == atom3 ) {
				my_ID = id2;
			} else if ( atom3->parent() == atom2 ) {
				my_ID = id3;
			} else  {
				utility_exit_with_message( "Problem with atoms: resno " + I(3,i) + " chino " + I(3,n) + " " +
																	 pose.residue(i).atom_name( id2.atomno() ) + " " +
																	 pose.residue(i).atom_name( id3.atomno() ) );
			}

			dof_id = DOF_ID( my_ID, D );
			std::cout << "Attempt to vary bond length for resno " << i << " atom: " << pose.residue(i).atom_name( my_ID.atomno() ) << std::endl;

			mm.set( dof_id, true );
			cst_set->add_dof_constraint( dof_id, new HarmonicFunc( (atom2->xyz() - atom3->xyz() ).length() , 0.1 ) );

		}
	}

	pose.constraint_set( cst_set );

}

/////////////////////////////////////////////////////////////////////////////////void
void
mask_out_loop( pose::Pose const & pose,
							 utility::vector1< std::map< Size, Size > > & alignment2sequence,
							 FArray1D_bool & sequence_mask )
{
	Size const alignment_length( alignment2sequence[1].size() );
	for (Size i = 1; i <= alignment_length; i++ ){
		if (  sequence_mask( i )  ) {
			Size const pdb_number( alignment2sequence[1][i] );
			if ( pose.secstruct( pdb_number ) == 'L' ) {
				sequence_mask( i ) = false;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
void vary_geometry_backbone( pose::Pose & pose, kinematics::MoveMap & mm )
{

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set( 	pose.constraint_set()->clone() ) ;

	//Change this to also include D DOF's for sidechains.
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		for (Size n = 1; n <= 3; n++ ) {
			TorsionID tor_id( i, BB, n );
			AtomID id1,id2,id3,id4, my_ID;
			bool const failure = pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );
			if (failure) continue;

			core::kinematics::tree::AtomCOP atom2 ( & pose.atom_tree().atom( id2 ) );
			core::kinematics::tree::AtomCOP atom3 ( & pose.atom_tree().atom( id3 ) );

			DOF_ID dof_id;
			if ( atom2->parent() == atom3 ) {
				my_ID = id2;
			} else if ( atom3->parent() == atom2 ) {
				my_ID = id3;
			} else  {
				utility_exit_with_message( "Problem with atoms: resno " + I(3,i) + " chino " + I(3,n) + " " +
																	 pose.residue(i).atom_name( id2.atomno() ) + " " +
																	 pose.residue(i).atom_name( id3.atomno() ) );
			}

			dof_id = DOF_ID( my_ID, D );
			std::cout << "Attempt to vary bond length for resno " << i << " atom: " << pose.residue(i).atom_name( my_ID.atomno() ) << std::endl;

			mm.set( dof_id, true );
			cst_set->add_dof_constraint( dof_id, new HarmonicFunc( (atom2->xyz() - atom3->xyz() ).length() , 0.1 ) );

		}

	// 	utility::vector1 < Size >  atom_indices( pose.residue(i).mainchain_atoms() );
// 		for (Size n = 1; n <= atom_indices.size(); n++ ) {

// 			DOF_ID dof_id;
// 			AtomID my_ID( atom_indices[ n ], i );

// 			dof_id = DOF_ID( my_ID, D );
// 			std::cout << "Attempt to vary bond length for resno " << i << " atom: " << pose.residue(i).atom_name( my_ID.atomno() ) << std::endl;

// 			mm.set( dof_id, true );
// 			pose.constraint_set()->add_dof_constraint( dof_id, new HarmonicFunc( pose.atom_tree().dof( dof_id ) , 0.2 ) );

// 		}
	}

	pose.constraint_set( cst_set );

}

/////////////////////////////////////////////////////////////////////////////////////////
void
setup_constraints( std::string const cst_file_name, pose::Pose & pose ){
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( cst_file_name, new ConstraintSet, pose );
// 	ConstraintSetOP cst_set_filtered = cst_set.clone();
// 	//Need to filter out constraints that don't apply:
// 	utility::vector1< ConstraintCOP>  const & all_constraints( cst_set.get_all_constraints() );
// 	for (Size n = 1; n <= all_constraints.size(); n++ ) {
// 		ConstraintCOP & cst( all_constraints[ n ] );
// 		for ( Size i=1; i<= cst->natoms(); ++i ) {
// 			int const seqpos( cst->atom(i).rsd() );

// 		}
// 	}

	pose.constraint_set( cst_set );

}


/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// WARNING! This only works for atom pairs!!!
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void
output_constraints( pose::Pose const & pose, std::string const & filename, std::map< Size, Size > & alignment2full ){
	using namespace core::scoring::constraints;
	utility::io::ozstream out( filename );
	//	ConstraintIO::get_instance()->write_constraints( out, pose );

	utility::vector1< ConstraintCOP>  const & all_constraints( pose.constraint_set()->get_all_constraints() );

	out << "[ atompairs ]" << std::endl;
 	for (Size n = 1; n <= all_constraints.size(); n++ ) {
 		ConstraintCOP const & cst( all_constraints[ n ] );

		if ( cst->natoms() != 2 ) continue;
 		for ( Size i=1; i<= cst->natoms(); ++i ) {
 			Size const seqpos( cst->atom(i).rsd() );
			out << A(5, pose.residue( seqpos ).atom_name( cst->atom(i).atomno() ) )
					<< ' ' << I(4, alignment2full[ seqpos ] );
 		}

		out << ' ';
		cst->get_func().show_definition( out );

 	}


}

////////////////////////////////////////////////////////////////////////////
void
output_constraints( pose::Pose const & pose, std::string const & filename ) {
	std::map< Size, Size > alignment2full;
	for (Size i=1; i <= pose.total_residue(); i++ ) alignment2full[ i ] = i;
	output_constraints( pose, filename, alignment2full );
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// This is lifted from Colin Smith's pilot app code (backrub.cc)
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void
backrub_protocol( pose::Pose const & native_pose, pose::PoseOP & pose, scoring::ScoreFunctionOP & score_fxn )
{

	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( new core::pack::task::operation::PreserveCBeta );

	score_fxn->set_weight(core::scoring::mm_bend, option[ backrub::mm_bend_weight ]);

	// set up the BackrubMover
	protocols::backrub::BackrubMover backrubmover;
	// read known and unknown optimization parameters from the database
	backrubmover.branchopt().read_database();

	// set up the SidechainMover
	protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover;
	sidechainmover.set_task_factory(main_task_factory);
	sidechainmover.set_prob_uniform(option[ backrub::sc_prob_uniform ]);

	//Wait a minute. SideChainMover is acting crazy. Need to set its task?
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( *pose ));
	task->initialize_from_command_line();
	for (Size i = 1; i <= pose->total_residue(); i++) {
	 	task->nonconst_residue_task(i).restrict_to_repacking();
	}

	sidechainmover.set_task( task );

	backrubmover.clear_segments();
	backrubmover.set_input_pose(pose);

	std::cout << "Backrub segment lengths: " << option[ backrub::min_atoms ] << "-" << option[ backrub::max_atoms ] << " atoms"
		 << std::endl;

	std::cout << "Backrub main chain pivot atoms: " << option[ backrub::pivot_atoms ].value_string() << std::endl;

	// determine list of residues to backrub
	utility::vector1<core::Size> resnums;
	if (option[ backrub::pivot_residues ].user()) {
		// if the user specified a vector of residues, get rid of any less than 1
		for (core::Size i = 1; i <= option[ backrub::pivot_residues ].size(); ++i) {
			if (option[ backrub::pivot_residues ][i] >= 1) resnums.push_back(option[ backrub::pivot_residues ][i]);
		}
	} else {
		// otherwise use all residues
		for (core::Size i = 1; i <= pose->total_residue(); ++i) resnums.push_back(i);
	}

	// add segments to the backrub mover
	backrubmover.add_mainchain_segments(resnums, option[ backrub::pivot_atoms ], option[ backrub::min_atoms ],
																			option[ backrub::max_atoms ]);

	std::cout << "Backrub Segments Added: " << backrubmover.num_segments() << std::endl;

	std::cout << "Score After PDB Load:" << std::endl;
	score_fxn->show(std::cout, *pose);

	backrubmover.optimize_branch_angles(*pose);
	sidechainmover.idealize_sidechains(*pose);

	std::cout << "Score After Branch Angle Optimization/Side Chain Idealization:" << std::endl;
	score_fxn->show(std::cout, *pose);

	Real const start_temperature = option[ backrub::mc_kt ];
	protocols::moves::MonteCarlo mc(*pose, *score_fxn, start_temperature);

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );


	using namespace core::scoring;

	kinematics::MoveMap mm;
	mm.set_jump( true );
	mm.set_bb( true );
	mm.set_chi( true );
	minimizer.run( *pose, mm, *score_fxn, options );
	score_fxn->show( std::cout, *pose );
	Real rmsd = CA_rmsd( native_pose, *pose );
	std::cout << "RMSD  ==> " << rmsd << std::endl;

	// create viewer windows if OpenGL is enabled
	//	protocols::viewer::add_monte_carlo_viewer(mc, "Backrub", 600, 600);

	// iterate to generate multiple structures
	// reset to the starting optimized pose
	//	core::pose::PoseOP pose(new core::pose::Pose(*pose));

	std::cout << "Running " << option[ backrub::ntrials ] << " trials..." << std::endl;

	static const Size NUM_ROUNDS( 20 );
	for (Size r = 1; r <= NUM_ROUNDS; r++ ) {

		mc.reset(*pose);
		mc.set_temperature( start_temperature * ( 1  - (r-1) * 1.0 /(NUM_ROUNDS-1) ) );

		for (int i = 1; i <= option[ backrub::ntrials ]; ++i) {

			std::string move_type;

			// could use random mover for this...
			if (RG.uniform() > option[ backrub::sc_prob ]) {
				backrubmover.apply(*pose);
				move_type = backrubmover.type();
			} else {
				sidechainmover.apply(*pose);
				move_type = sidechainmover.type();
			}

			mc.boltzmann(*pose, move_type);

			rmsd = CA_rmsd( native_pose, *pose );
			std::cout << "RMSD [" << I(4, i) << "] ==> " << rmsd << std::endl;
		}

		mc.show_counters();

		// dump out the low score and last accepted poses
		std::cout << "Last Score:" << std::endl;
		score_fxn->show(std::cout, *pose);

		std::cout << "Low Score:" << std::endl;
		score_fxn->show(std::cout, *pose);

		minimizer.run( *pose, mm, *score_fxn, options );
		score_fxn->show( std::cout, *pose );
		rmsd = CA_rmsd( native_pose, *pose );
	}

	*pose = mc.lowest_score_pose();

	pose->dump_pdb( "backrub_final.pdb" );

}

///////////////////////////////////////////////////////////////////////////////////////
void
minimize_protocol( pose::Pose const & native_pose, pose::Pose & pose, scoring::ScoreFunctionOP & scorefxn )
{

	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	Real init_rmsd = CA_rmsd( native_pose, pose );

	/////////////////////////////////////////////////////////////////////////////
	// Minimize
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	Real rmsd;
	mm.set_jump( false );
	mm.set_bb( false );
	mm.set_chi( false );

	static bool const VARY_GEOMETRY( option[ vary_geometry]() );

	mm.set_chi( true );
	if ( VARY_GEOMETRY ) vary_geometry_sidechains( pose, mm );
	scorefxn->set_weight( dof_constraint, 1.0 );
	scorefxn->show( std::cout, pose );

	minimizer.run( pose, mm, *scorefxn, options );
	scorefxn->show( std::cout, pose );
	rmsd = CA_rmsd( native_pose, pose );
	std::cout << "RMSD " << init_rmsd << " [init]   " << rmsd << " [final] " << std::endl;
	pose.dump_pdb( "minimize_chi.pdb" );

	mm.set_jump( true );
 	minimizer.run( pose, mm, *scorefxn, options );
 	scorefxn->show( std::cout, pose );
 	rmsd = CA_rmsd( native_pose, pose );
 	std::cout << "RMSD " << init_rmsd << " [init]   " << rmsd << " [final] " << std::endl;
 	pose.dump_pdb( "minimize_chi_jump.pdb" );

	mm.set_bb(  true );
	if ( VARY_GEOMETRY ) 	vary_geometry_backbone( pose, mm );
	minimizer.run( pose, mm, *scorefxn, options );
	scorefxn->show( std::cout, pose );
	rmsd = CA_rmsd( native_pose, pose );
	std::cout << "RMSD " << init_rmsd << " [init]   " << rmsd << " [final] " << std::endl;
	pose.dump_pdb( "minimize_chi_jump_bb.pdb" );

	std::cout << "*************************************************" << std::endl;
	std::cout << "DONE WITH QUICK MIN ==> minimize_chi_jump_bb.pdb " << std::endl;
	std::cout << "*************************************************" << std::endl;

	Size const num_rounds( option[ rounds ]() );

	for (Size n = 1; n <= num_rounds; n++ ) {
		repack( pose, *scorefxn );
		minimizer.run( pose, mm, *scorefxn, options );
		scorefxn->show( std::cout, pose );
		rmsd = CA_rmsd( native_pose, pose );
		std::cout << "RMSD " << init_rmsd << " [init]   " << rmsd << " [final] " << std::endl;
		pose.dump_pdb( "minimize_chi_jump_bb_afterpack.pdb" );
	}

	if (false){
		std::cout << "Turning off CA tether " << std::endl;
		scorefxn->set_weight( atom_pair_constraint, 0.0  );
		minimizer.run( pose, mm, *scorefxn, options );
		scorefxn->show( std::cout, pose );
		rmsd = CA_rmsd( native_pose, pose );
		std::cout << "RMSD " << init_rmsd << " [init]   " << rmsd << " [final] " << std::endl;
		pose.dump_pdb( "minimize_chi_jump_bb_notether.pdb" );
	}


}


////////////////////////////////////////////////////////////////////////////
void
output_constraints_for_full_length( pose::Pose const & pose,
																		std::map< Size, Size > & alignment2sequence_target,
																		FArray1D_bool const & sequence_mask,
																		std::string const & filename)
{
	//Need to figure out mapping from this truncated pose sequence to full-length sequence
	std::map< Size, Size > alignment2full;
	Size count( 0 );
	for (Size i = 1; i <= alignment2sequence_target.size(); i++ ){
		if (  sequence_mask( i )  ) {
			count++;
			Size const pdb_number( alignment2sequence_target[i] );
			alignment2full[ count ] = pdb_number;
		}
	}
	output_constraints( pose, filename, alignment2full );
}

////////////////////////////////////////////////////////////////////////////
void
easy_target_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;

	////////////////////////////////////////////////
	//Read in fasta file.
	std::string const fasta_file ( option[ in::file::fasta ]() );
	utility::vector1< std::string > sequences, pdb_names;
	Size const alignment_length = read_alignment_fasta_file( sequences, pdb_names, fasta_file );

	//Figure out a mask.
	FArray1D_bool sequence_mask( alignment_length, false );
	setup_mask( sequence_mask, sequences );

	//Get from alignment to index inside each pdb.
	utility::vector1< std::map< Size, Size > > alignment2sequence;
	setup_alignment_map( alignment2sequence, sequences );

	////////////////////////////////////////////////
	//Read in native pdb
	Pose native_pose, template_pose, temp_pose;
	PoseOP pose_op( new Pose );
	Pose & pose( *pose_op );

	std::string native_file = option[ in::file::native ];
	io::pdb::pose_from_pdb( native_pose, native_file );

	//Read in template pdb.
	std::string template_file = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( template_pose, template_file );

	//Read in secstruct?
	//	std::string template_secstruct_file = option[ secstruct_file ];
	//	setup_secstruct( template_pose, template_secstruct_file );
	setup_secstruct_dssp( template_pose );

	if ( option[ mask_loop ] ){
		mask_out_loop( template_pose, alignment2sequence, sequence_mask );
	}

	// Apply mask to native and to template.
	std::string masked_sequence = apply_mask( sequence_mask, alignment2sequence, pdb_names, native_file, native_pose );
	apply_mask( sequence_mask, alignment2sequence, pdb_names, template_file, template_pose );

	native_pose.dump_pdb( "native.pdb" );
	template_pose.dump_pdb( "template.pdb" );


	//	Pose pose_with_desired_sequence = native_pose; //For now.
	//	Size const nres_model = pose_with_desired_sequence.total_residue();

	std::string desired_sequence = "";
	for (Size n = 1; n <= sequences[1].length(); n++ ) {
		if( sequences[1][n-1] != '-' ) desired_sequence += sequences[1][n-1];
	}
	Size const nres_model = desired_sequence.size();

	std::cout << "WANT SEQUENCE! ==> " << std::endl;
	std::cout << desired_sequence << std::endl;

	//Prepare template pdb.
	FArray1D_bool conserved( nres_model, false );
	FArray1D_bool sequence_mask_local( nres_model, true );

	utility::io::ozstream out( "chi_stats.txt" );

	ScoreFunctionOP scorefxn = getScoreFunction();
	(*scorefxn)( native_pose ); //for neighbor info. This seems silly.

	//	prepare_start_model( template_pose, pose_with_desired_sequence, *scorefxn, pose, conserved );
	prepare_start_model( template_pose, desired_sequence, *scorefxn, pose, conserved );

	pose.dump_pdb( "start.pdb" );
	output_chi_stats( native_pose, pose, sequence_mask_local, out, "standard.wts", 1 );

	output_backbone_stats( native_pose, pose );

	/////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << "Conserved " << std::endl;
	//		output_chi_stats( native_pose, template_pose, conserved, out, "conserved template", 2 );
	//		output_chi_stats( native_pose, pose, conserved, out, "conserved repack", 3 );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//Constraints?
	if ( option[ cst_file ].user() ) {
		setup_constraints( option[ cst_file ], pose );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}

	if ( option[ CA_tether ].user() ) {
		if ( option[ use_native_CA ] ) {
			setup_CA_constraints( pose, native_pose );
		}	else {
			setup_CA_constraints( pose, pose);
		}
		scorefxn->set_weight( atom_pair_constraint, 1.0  );
	}

	if ( option[ rhiju_fold_tree ] ) {
		setup_rhiju_fold_tree( pose );
		scorefxn->set_weight( chainbreak, 1.0 );
	}

	output_constraints( pose, "new.constraints" );

	output_constraints_for_full_length( pose, alignment2sequence[1], sequence_mask, "new_full_length.constraints" );

	static bool const try_backrub = option[ backrub_test ];

	if ( try_backrub ) {
		backrub_protocol ( native_pose, pose_op, scorefxn );
	} else {
		minimize_protocol( native_pose, pose, scorefxn );
	}


// 	if (false)
// 	{
// 		ScoreFunctionOP scorefxn = getScoreFunction();
// 		scorefxn->set_weight( fa_rep, 0.1 );
// 		prepare_start_model( template_pose, pose_with_desired_sequence, *scorefxn, pose, conserved );
// 		pose.dump_pdb( "start_low_farep.pdb" );
// 		output_chi_stats( native_pose, pose, sequence_mask_local, out, "low_fa_rep", 2 );
// 	}

// 	if (false)
// 	{
// 		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( SOFT_REP_WTS );
// 		prepare_start_model( template_pose, pose_with_desired_sequence, *scorefxn, pose, conserved );
// 		pose.dump_pdb( "start_soft_rep.pdb" );
// 		output_chi_stats( native_pose, pose, sequence_mask_local, out, "soft_rep.wts", 3 );
// 	}

}

void
prepare_full_length_start_model(
	 pose::Pose & template_pose,
	 pose::Pose & pose,
	 utility::vector1<std::string> const & sequences,
	 FArray1D_bool const & sequence_mask,
	 utility::vector1< std::map< Size, Size > > & alignment2sequence,
	 utility::vector1< std::string > const & pdb_names,
	 std::string const & which_file
																)
{
	using namespace core::chemical;
	using namespace core::conformation;

	pose.clear();
	Size const alignment_length( alignment2sequence[1].size() );

	// Find which alignment to use, based on input pdb name.
	Size which_sequence( 0 );
	for( Size n=1; n<=pdb_names.size(); n++ ){
		//		std::cout << pdb_names[n] << " " << which_file << std::endl;
		if ( pdb_names[n] == which_file ) {
			which_sequence = n;
			break;
		}
	}

	ResidueTypeSet const & rsd_set( template_pose.residue(1).residue_type_set() );

	std::string const full_desired_sequence = sequences[1];

	std::string desired_sequence = "";
	for (Size i = 1; i <= alignment_length; i++ ){
		if ( !sequence_mask( i ) )  continue;
		Size const pdb_number( alignment2sequence[ which_sequence ][i] );
		pose.append_residue_by_bond( template_pose.residue( pdb_number ) );
		desired_sequence += desired_sequence[i-1];
	}

	//Write over sequence?
	for (Size i = 1; i <= desired_sequence.size(); i++ ){
		char const new_seq = desired_sequence[i-1];
		ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );
		pose.replace_residue( i, *new_rsd, false );
	}

	//Now need to add in residues that don't exist already. This is tricky. Need to do it segment by segment.
	//Loop definition
	bool in_loop( false );
	Size start( 0 ), end( 0 );
	utility::vector1< std::pair< Size, Size > > loops;
	for (Size i = 1; i <= alignment_length; i++ ){
		Size const pdb_number( alignment2sequence[ 1 ][i] );
		if (pdb_number == 0 ) continue;
		//		std::cout << i <<  " " << pdb_number << " " << sequence_mask(i) << " " << in_loop << std::endl;

		if ( !sequence_mask(i) && !in_loop) {
			in_loop = true;
			start = pdb_number;
		}
		if ( sequence_mask(i) && in_loop ) {
			in_loop = false;
			end = pdb_number-1;
			loops.push_back( std::make_pair( start, end ) );
		}
	}


	//1. Will need to handle special case for termini -- not coded yet!!
	//2. Probably will want to randomly move around cutpoint (by "prepending"!).
	for (Size n = 1; n <= loops.size(); n++ ) {
		std::cout << "LOOP " << n << " ==> " << loops[n].first << " " << loops[n].second << std::endl;
		for (Size i = loops[n].first; i <= loops[n].second; i++ ) {
			ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( desired_sequence[i-1] ).exclude_variants().select( rsd_set )[1] );
			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
			pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, i-1, true );
			pose.set_omega( i, 180.0 );
		}
	}

}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void
easy_loop_model_test(){

	//Gee whiz, I think I might need to write this from scratch!
	// Then worry about integrating with the protocol above.
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;

	////////////////////////////////////////////////
	//Read in fasta file.
	std::string const fasta_file ( option[ in::file::fasta ]() );
	utility::vector1< std::string > sequences, pdb_names;
	Size const alignment_length = read_alignment_fasta_file( sequences, pdb_names, fasta_file );

	//Figure out a mask.
	FArray1D_bool sequence_mask( alignment_length, false );
	setup_mask( sequence_mask, sequences );

	//Get from alignment to index inside each pdb.
	utility::vector1< std::map< Size, Size > > alignment2sequence;
	setup_alignment_map( alignment2sequence, sequences );

	////////////////////////////////////////////////
	//Read in native pdb
	Pose native_pose, template_pose, temp_pose;
	PoseOP pose_op( new Pose );
	Pose & pose( *pose_op );

	std::string native_file = option[ in::file::native ];
	io::pdb::pose_from_pdb( native_pose, native_file );

	//Read in template pdb.
	std::string template_file = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( template_pose, template_file );

	prepare_full_length_start_model( template_pose, pose, sequences, sequence_mask,  alignment2sequence, pdb_names, template_file );

	pose.dump_pdb( "start.pdb" );

	//////////////////////////////////////////////////////////////
	// Prepare starting pose. with extended, ideal loops. Hmmmm.
	// Copy from phil's test1.cc?

	//////////////////////////////////////////////////////////////

	//	ScoreFunctionOP scorefxn = getScoreFunction();
	//	protocols::LoopRebuild loop_builder( *scorefxn, loops, frag_libs );//Will this work as fullatom?
	//	loop_rebuilder.apply( pose );

	//	protocols::LoopRemodelMover loopremodel_mover;



}


////////////////////////////////////////////////////////////////////////////////////////
void
initialize_sequence_mask( pose::Pose & pose, FArray1D_bool & sequence_mask ) {

	using namespace core::options;
	using namespace core::options::OptionKeys;

	sequence_mask.dimension( pose.total_residue() );
	sequence_mask = true;

	if ( option[ sequence_mask_file ].user() ) {
		utility::io::izstream data_stream( option[ sequence_mask_file ] );
		std::string line;
		getline(data_stream, line);
		for (Size n = 1; n <= pose.total_residue(); n++ ) {
			sequence_mask( n ) = line[n-1];
		}
	}

}
////////////////////////////////////////////////////////////////////////////////////////
void
chi_stats_test()
{
	using namespace core::options;
	using namespace core::scoring;

	pose::Pose native_pose;
	std::string native_file = option[ in::file::native ];
	io::pdb::pose_from_pdb( native_pose, native_file );

	///////////////////////////////////
	// To setup neighbors, burial info.
	ScoreFunctionOP scorefxn = getScoreFunction()  );
	( *scorefxn )( native_pose );

	///////////////////////////////////////////////////
	// In case we don't want to look at all residues.
	FArray1D_bool sequence_mask;
	initialize_sequence_mask( native_pose, sequence_mask );

	utility::vector1< std::string > pdb_files = option[ in::file::s ]();

	std::string outfile = option [ out::file::o ]();
	utility::io::ozstream out( outfile );

	///////////////////////////////////
	for (Size n = 1; n <= pdb_files.size(); n++ ) {

		pose::Pose pose;
		std::string const pdb_file = pdb_files[n];

		std::cout << " About to read in: " << pdb_file << std::endl;
		io::pdb::pose_from_pdb( pose, pdb_file );

		//if ( pose.sequence() != native_pose.sequence() ) continue;
		output_chi_stats( native_pose, pose, sequence_mask, out, pdb_file, n );
	}


}

//////////////////////////////////////////////////////////////////////////////////////////////////
void
cst_relax_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;

	utility::vector1< std::string > start_files ( option[ in::file::s ]() );

	PoseOP native_pose_op( new pose::Pose ) ;
	bool use_native( false );
	if ( option[ in::file::native  ].user() ) {
		std::string native_file = option[ in::file::native ];
		io::pdb::pose_from_pdb( *native_pose_op, native_file );
		use_native = true;
	}


	ScoreFunctionOP scorefxn = getScoreFunction()  );

	// Basic set up
	bool const score_only = option[ run::score_only ]();
	protocols::simple_filters::RmsdEvaluator rmsd_evaluator( native_pose_op );
	rmsd_evaluator.report_gdt_components( true );
	core::io::silent::SilentFileData silent_file_data;
	BinarySilentStruct s;
	Size const nstruct = option[ out::nstruct ]();
	std::string const out_path = option[ out::path::path ]();
	std::string const silent_file = out_path+'/'+option[ out::file::silent  ]();


	//READ IN FIRST POSE.
	Pose pose, input_pose;
	io::pdb::pose_from_pdb( input_pose, start_files[1] );
	pose = input_pose;
	//graphics!
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	for (Size j = 1; j <= start_files.size(); j++ ) {

		if ( j > 1 ) io::pdb::pose_from_pdb( input_pose, start_files[j] );
		//		pose.dump_pdb( "start.pdb" );

		// Main loop
		for (Size n = 1; n <= nstruct; n++ ) {

			pose = input_pose;
			if ( option[ cst_file ].user() ) {
				setup_constraints( option[ cst_file ], pose );
				scorefxn->set_weight( atom_pair_constraint, 1.0 );
			}

			(*scorefxn)(pose); //score it.

			if ( ! score_only ) {
				relax::ClassicRelax relax_protocol( scorefxn );
				//	relax_protocol.set_stage2_cycles( 25 );
				//	relax_protocol.set_stage2_repack_period( 10 );
				//	relax_protocol.set_stage3_cycles( 10 );
  			relax_protocol.apply( pose );
			}

			//		setPoseExtraScores( pose, "rms",   protocols::simple_filters::native_CA_rmsd(native_pose, pose ) );

			std::string out_file_tag =   "S_"+lead_zero_string_of( n, 4 );
			if ( score_only ) out_file_tag = start_files[j];

			std::cout << "Making silent struct for " << out_file_tag << std::endl;
			s.fill_struct( pose,  out_file_tag );//, true /*fullatom*/ );
			if ( use_native ) rmsd_evaluator.apply( pose, out_file_tag, s );
			silent_file_data.write_silent_struct( s, silent_file );//, true /*write score only*/ );
		//		pose.dump_pdb( out_path + '/' + out_file_tag+".pdb" );

		}

	}


}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ chi_stats]() ) {
		chi_stats_test();
	} else if ( option[ easy_loop_model ]() ) {
		easy_loop_model_test();
	} else if ( option[ cst_relax ]() ) {
		cst_relax_test();
	} else {
		easy_target_test();
	}

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( chi_stats, "", false );
	NEW_OPT( copy_native_chi, "", false );
	NEW_OPT( rhiju_fold_tree, "", false );
	NEW_OPT( use_native_CA, "", false );
	NEW_OPT( vary_geometry, "", false );
	NEW_OPT( mask_loop, "", false );
	NEW_OPT( backrub_test, "", false );
	NEW_OPT( easy_loop_model, "", false );
	NEW_OPT( cst_relax, "", false );
	NEW_OPT( t469_zinc_tether, "", false );
	NEW_OPT( soft_segment_file, "", "blah.txt" );
	NEW_OPT( sequence_mask_file, "", "sequence_mask.txt" );
	NEW_OPT( secstruct_file, "", "blah.secstruct" );
	NEW_OPT( cst_file, "", "cst.txt" );
	NEW_OPT( rounds, "", 1 );
	NEW_OPT( cst_trim_loop, "", 2 );
	NEW_OPT( chi1_constraint_weight, "", 0.0 );
	NEW_OPT( chi2_constraint_weight, "", 0.0 );
	NEW_OPT( CA_tether, "", 0.0 );
	NEW_OPT( soft_CA_tether, "", 5.0 );

	//From colin's backrub.cc
	NEW_OPT(backrub::pivot_residues, "residues for which contiguous stretches can contain segments (internal residue numbers, defaults to all residues)", utility::vector1<int>());
	NEW_OPT(backrub::pivot_atoms, "main chain atoms usable as pivots", utility::vector1<std::string>(1, "CA"));
	NEW_OPT(backrub::min_atoms, "minimum backrub segment size (atoms)", 3);
	NEW_OPT(backrub::max_atoms, "maximum backrub segment size (atoms)", 34);
	NEW_OPT(backrub::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(backrub::sc_prob, "probability of making a side chain move", 0.25);
	NEW_OPT(backrub::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.1);
	NEW_OPT(backrub::mc_kt, "value of kT for Monte Carlo", 0.3);
	NEW_OPT(backrub::mm_bend_weight, "weight of mm_bend bond angle energy term", 1.0);

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
