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
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/scoring/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <protocols/farna/util.hh>

#include <protocols/viewer/viewers.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//StepWise!
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/protein/util.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pose_stream/ExtendedPoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <core/util/datacache/BasicDataCache.hh>
#include <core/util/datacache/CacheableString.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>


#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
//RNA stuff.
//#include <protocols/farna/fragments/RNA_FragmentsClasses.hh>
//#include <protocols/farna/RNA_DeNovoProtocol.hh>
//#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <vector>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, cluster_test )
OPT_KEY( StringVector, s1 )
OPT_KEY( StringVector, s2 )
OPT_KEY( StringVector, silent1 )
OPT_KEY( StringVector, silent2 )
OPT_KEY( StringVector, tags1 )
OPT_KEY( StringVector, tags2 )
OPT_KEY( IntegerVector, slice_res1 )
OPT_KEY( IntegerVector, slice_res2 )
OPT_KEY( IntegerVector, input_res1 )
OPT_KEY( IntegerVector, input_res2 )
OPT_KEY( String, cst_file )
OPT_KEY( String, pack_weights )
OPT_KEY( String, align_pdb )
OPT_KEY( Boolean, cluster_by_all_atom_rmsd )
OPT_KEY( Boolean, output_start )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, parse_pathway )


// ///////////////////////////////////////////////////////////////////////////////
// // might be useful in a util.cc somewhere
// core::scoring::constraints::ConstraintSetOP
// constraint_set_slice( core::scoring::constraints::ConstraintSetOP & cst_set,
// 												utility::vector1< core::Size > const & slice_res )
// {

// 	using namespace core::scoring::constraints;
// 	using namespace core::scoring;
// 	using namespace core::id;

// 	ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

// 	ConstraintCOPs csts( cst_set->get_all_constraints() );

// 	std::map< Size, Size > slice_map;
// 	for (Size i = 1; i <= slice_res.size();i++) slice_map[ slice_res[ i ] ] = i;

// 	for ( Size n = 1; n <= csts.size(); n++ ) {

// 		ConstraintCOP const & cst( csts[n] );

// 		if ( cst->score_type() == atom_pair_constraint)  { // currently only defined for pairwise distance constraints.
// 			Size const i = cst->atom( 1 ).rsd();
// 			Size const j = cst->atom( 2 ).rsd();
// 			//			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
// 			//			if ( dist  > separation_cutoff ) continue;

// 			if ( slice_map.find( i ) == slice_map.end()  ) continue;
// 			if ( slice_map.find( j ) == slice_map.end()  ) continue;

// 			//std::cout << "CST MAP: " << i << " " << slice_map[ i] << "          " << j << " " << slice_map[ j ] << std::endl;
// 			AtomID atom1_new( cst->atom(1).atomno(),  slice_map[ i ]);
// 			AtomID atom2_new( cst->atom(2).atomno(),  slice_map[ j ]);
// 			ConstraintOP cst_new = new AtomPairConstraint( atom1_new, atom2_new,
// 																										 cst->get_func().clone() /*is this defined?*/, cst->score_type() );

// 			cst_set_new->add_constraint( cst_new );

// 		}

// 	}


// 	std::cout << "NUM CONSTRAINTS " << cst_set_new->get_all_constraints().size() << " out of " <<
// 		csts.size() << std::endl;

// 	return cst_set_new;
// }


///////////////////////////////////////////////////////////////////////////////
core::io::pose_stream::PoseInputStreamOP
setup_pose_input_stream(
												utility::vector1< std::string > const & option_s1,
												utility::vector1< std::string > const & option_silent1,
												utility::vector1< std::string > const & option_tags1
												){
	using namespace core::io::pose_stream;

	PoseInputStreamOP input1;

	if( option_s1.size() > 0 ) {
		// pdb input(s).
		input1 = new PDBPoseInputStream( option_s1 );

	} else if ( option_silent1.size() > 0 ){

		if ( option_tags1.size() > 0) {
			input1 = new SilentFilePoseInputStream( option_silent1 ,
																							option_tags1 );
		} else {
			input1 = new SilentFilePoseInputStream( option_silent1 );
		}
	} else {
		// create a pose stream with a single blank pose...
		input1 = new ExtendedPoseInputStream( "", 1 ); // hmm...
	}

	return input1;

}

///////////////////////////////////////////////////////////////
void
parse_pathway_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace core::pose;
	using namespace core::pack;
	using namespace protocols::swa;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	Pose start_pose;
	std::string const input_pdb_file = option[ in::file::s]()[1];
	io::pdb::pose_from_file( start_pose, *rsd_set, input_pdb_file , core::import_pose::PDB_file);

	Size const nres( start_pose.total_residue() );

	ScoreFunctionOP minimize_scorefxn = get_score_function();

	// initialize minimizer
	bool const deriv_check_( false /*option[ deriv_check]*/ );
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	mm.set_jump( true );

	std::cout << "Initialize minimization... " << std::endl;
	minimizer.run( start_pose, mm, *minimize_scorefxn, options );
	minimize_scorefxn->show( std::cout, start_pose );

	std::string silent_file = input_pdb_file + ".sc";
	if( option[ out::file::silent  ].user() )  silent_file = option[ out::file::silent  ]();

	bool const add_peptide_plane_( true );

	for ( Size i = 1; i <= nres; i++ ) {
		for ( Size j = i; j <= nres; j++ ) {

			Pose pose = start_pose;

			utility::vector1 < Size > slice_res;
			for ( Size n = i; n <= j; n++ ) slice_res.push_back( n );
			protocols::stepwise::pdbslice( pose, slice_res );

			if ( i > 1 && add_peptide_plane_ ) {
				//Pose pose_copy = pose;
				//remove_variant_type_from_pose_residue( pose, "LOWER_TERMINUS", 1  );
				//				pose.set_xyz( core::id::AtomID( pose.residue(1).atom_index("H"),  1 ),
				//											pose_copy.xyz( core::id::AtomID( pose.residue( 1 ).atom_index( "1H" ), 1 ) ) );
				add_variant_type_to_pose_residue( pose, "N_ACETYLATION", 1 );
				pose.set_phi( 1, start_pose.phi( i ) ) ;
			}


			if ( j < nres && add_peptide_plane_ ) {
				//				remove_variant_type_from_pose_residue( pose, "UPPER_TERMINUS", pose.total_residue()  );
				add_variant_type_to_pose_residue( pose, "C_METHYLAMIDATION", pose.total_residue() );
				pose.set_psi( pose.total_residue(), start_pose.psi( j ) ) ;
				pose.set_omega( pose.total_residue(), start_pose.omega( j ) ) ;
			}

			std::cout << i << " " << j << std::endl;
			(*minimize_scorefxn)( pose );

			std::string const tag = "S_"+string_of(i)+"_"+string_of(j);
			BinarySilentStruct s( pose, tag );
			s.add_energy( "start", i );
			s.add_energy( "end", j );

			static const SilentFileData silent_file_data;
			silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );

		}
	}


}

///////////////////////////////////////////////////////////////////////////////
void
stepwise_template_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace core::pose;
	using namespace core::pack;
	using namespace protocols::stepwise::protein;
	using namespace protocols::swa;

	//////////////////////////
	// read in total sequence
	//Read in desired fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();


	//////////////////////
	// POSE SETUP
	//////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Pose full_pose;
	make_pose_from_sequence( full_pose, desired_sequence, *rsd_set );
	bool read_start_pdb_from_disk( false );
	//	if ( option[ in::file::s ].user() ) {
		//		io::pdb::pose_from_file( full_pose, *rsd_set, option[ in::file::s]()[1] , core::import_pose::PDB_file);
		//		read_start_pdb_from_disk = true;
	//	}

	/////////////////////
	// Read in native
	PoseOP native_pose, align_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_file( *native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
		if ( desired_sequence != native_pose->sequence() ) utility_exit_with_message( "Native pose sequence looks wrong." );
	}
	if (option[ align_pdb ].user() ) {
		align_pose = PoseOP( new Pose );
		std::string align_pdb_file  = option[ align_pdb ]();
		io::pdb::pose_from_file( *align_pose, *rsd_set, align_pdb_file , core::import_pose::PDB_file);
	}


	///////////////////////
	// What's coming in?
	utility::vector1< Size > const & input_res1_ = option[ input_res1 ]();
	utility::vector1< Size > const & input_res2_ = option[ input_res2 ]();

	///////////////////////////////////////////
	// what sequence are we working with?
	Size const totres = desired_sequence.size();
	ObjexxFCL::FArray1D< bool > is_input_res( totres, false );
	std::map < Size, Size > full_to_sub, res_map1, res_map2;
	for ( Size i = 1; i <= input_res1_.size(); i++ ) is_input_res( input_res1_[ i ] ) = true;
	for ( Size i = 1; i <= input_res2_.size(); i++ ) is_input_res( input_res2_[ i ] ) = true;

	std::string desired_sub_sequence = "";
	utility::vector1< Size > working_res, superimpose_res;
	Size count( 0 );
	for ( Size i = 1; i <= totres; i++ ) {
		if ( is_input_res( i ) ) {
			count += 1;
			desired_sub_sequence += desired_sequence[i-1];
			working_res.push_back( i );
			full_to_sub[ i ] = count;
			superimpose_res.push_back( count );
		}
	}

	for ( Size i = 1; i <= input_res1_.size(); i++ )  res_map1[ full_to_sub[ input_res1_[ i ] ] ] = i;
	for ( Size i = 1; i <= input_res2_.size(); i++ )  res_map2[ full_to_sub[ input_res2_[ i ] ] ] = i;

	///////////////////////////////////////////
	// initial pose setup.
	Pose pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	std::cout << "DESIRED SUB SEQUENCE " << desired_sub_sequence << std::endl;
	pose = full_pose;
	pdbslice( pose, working_res );
	if ( native_pose ) pdbslice( *native_pose, working_res );
	if ( align_pose )  pdbslice( *align_pose,  working_res );

	// make extended chain
	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		if ( ! pose.residue(pos).is_protein() ) continue;
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );
	}
	//	pose.dump_pdb( "extended.pdb" );

	// Constraints...
	ConstraintSetOP cst_set;
	bool constraints_exist( false );
	if ( option[ cst_file ].user() ) {
		cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, full_pose );
		cst_set = constraint_set_slice( cst_set, working_res, pose, full_pose );
		pose.constraint_set( cst_set );
		constraints_exist = true;
	}

	//////////////////////////////////////////////////////
	//  packer and minimizer setup.
	ScoreFunctionOP pack_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
	ScoreFunctionOP minimize_scorefxn = get_score_function();

	if ( constraints_exist ){
		pack_scorefxn->set_weight( atom_pair_constraint, 1.0 );
		minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}

	// initialize packer
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->restrict_to_repacking();

	// initialize minimizer
	bool const deriv_check_( false /*option[ deriv_check]*/ );
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );
	mm.set_jump( true );


	//////////////////////////////////////////
 	// main loop.
	//////////////////////////////////////////
	PoseInputStreamOP stream1 = setup_pose_input_stream( option[ s1 ](), option[ silent1 ](), option[ tags1 ]() );
	PoseInputStreamOP stream2 = setup_pose_input_stream( option[ s2 ](), option[ silent2 ](), option[ tags2 ]() );


	std::string const silent_file = option[ out::file::silent  ]();

	utility::vector1< Size > const slice_res1_ = option[ slice_res1 ]();
	utility::vector1< Size > const slice_res2_ = option[ slice_res2 ]();
 	Pose pose_silent, input_pose1, input_pose2;
	Pose start_pose = pose;
	Size count1( 0 ), count2( 0 );
 	while ( stream1->has_another_pose() ) {

		pose = start_pose;

		stream1->fill_pose( input_pose1, *rsd_set );
		//		pose.dump_pdb( "start0.pdb" );

		if ( input_pose1.total_residue() > 0 ) {
			if ( slice_res1_.size() > 0 ) pdbslice( input_pose1, slice_res1_ );
			//input_pose1.dump_pdb( "input_pose1.pdb" );
			std::cout << pose.annotated_sequence( true ) << std::endl;
			std::cout << input_pose1.annotated_sequence( true ) << std::endl;
			copy_dofs( pose, input_pose1, res_map1 );
		}
		count1 += 1;

		//		pose.dump_pdb( "start1.pdb" );

		while ( stream2->has_another_pose() ) {

			stream2->fill_pose( input_pose2, *rsd_set );
			if ( input_pose2.total_residue() > 0 ) {
				if ( slice_res2_.size() > 0 ) 	 pdbslice( input_pose2, slice_res2_ );
				copy_dofs( pose, input_pose2, res_map2 );
			}
			count2 += 1;

			//			pose.dump_pdb( "start2.pdb" );

			(*minimize_scorefxn)( pose );
			if ( option[ output_start ] ) output_silent_struct( pose, native_pose, silent_file, "START" );

			if ( align_pose ){
				id::AtomID_Map< id::AtomID > const & alignment_atom_id_map =
					create_alignment_id_map( pose, *align_pose, superimpose_res );
				core::scoring::superimpose_pose( pose, *align_pose, alignment_atom_id_map);
			}


			pack_scorefxn->show( std::cout, pose );
			pack::pack_rotamers( pose, *pack_scorefxn, task);
			pack_scorefxn->show( std::cout, pose );
			minimize_scorefxn->show( std::cout, pose );
			minimizer.run( pose, mm, *minimize_scorefxn, options );
			minimize_scorefxn->show( std::cout, pose );

			//     save to silent file.
			output_silent_struct( pose, native_pose, silent_file, "S_"+lead_zero_string_of( count1, 4 ) + "_" + lead_zero_string_of( count2, 4) );

		}

	}

}


///////////////////////////////////////////////////////////////
void
cluster_outfile_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	protocols::stepwise::StepWiseLegacyClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	stepwise_clusterer.set_cluster_radius(	option[ OptionKeys::cluster::radius ]()	);
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] );
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );
	stepwise_clusterer.set_auto_tune( option[ auto_tune ] );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ cluster_test ] ){
		cluster_outfile_test();
	} else if ( option[ parse_pathway ] ){
		parse_pathway_test();
	} else {
		stepwise_template_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	//Uh, options?
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );

	NEW_OPT( s1, "input file(s)", blank_string_vector );
	NEW_OPT( s2, "input file(s)", blank_string_vector );
	NEW_OPT( silent1, "input file", blank_string_vector );
	NEW_OPT( silent2, "input file", blank_string_vector );
	NEW_OPT( tags1, "input tag(s)", blank_string_vector );
	NEW_OPT( tags2, "input tag(s)", blank_string_vector );
	NEW_OPT( slice_res1, "Residues to slice out of starting file", blank_size_vector );
	NEW_OPT( slice_res2, "Residues to slice out of starting file", blank_size_vector );
	NEW_OPT( input_res1, "Residues already present in starting file", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting file2", blank_size_vector );
	NEW_OPT( pack_weights, "weights for green packing", "standard.wts" );
	NEW_OPT( align_pdb, "PDB file for alignment", "" );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( cluster_by_all_atom_rmsd, "cluster by all atom rmsd", false );
	NEW_OPT( output_start, "output starting pdb", false );
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( parse_pathway, "parse the pathway", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
