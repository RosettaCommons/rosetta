// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief


// libRosetta headers
 #include <protocols/frags/VallData.hh>
 #include <protocols/frags/TorsionFragment.hh>

 #include <core/scoring/dna/setup.hh>
 #include <core/scoring/dna/base_geometry.hh>
 #include <core/scoring/dna/BasePartner.hh>
 #include <core/scoring/GenBornPotential.hh>
 #include <core/scoring/LREnergyContainer.hh>
 #include <core/scoring/methods/Methods.hh>

 #include <protocols/simple_moves/BackboneMover.hh>
 #include <protocols/simple_moves/MinMover.hh>
 #include <protocols/moves/MonteCarlo.hh>
 #include <protocols/moves/Mover.hh>
 #include <protocols/moves/MoverContainer.hh>
 #include <protocols/moves/OutputMovers.hh>
 #include <protocols/rigid/RigidBodyMover.hh>
 // #include <protocols/moves/rigid_body_moves.hh>
 #include <protocols/moves/TrialMover.hh>
 #include <protocols/simple_moves/PackRotamersMover.hh>
 #include <protocols/simple_moves/RotamerTrialsMover.hh>
 #include <protocols/moves/RepeatMover.hh>

 #include <protocols/loops/ccd_closure.hh>
 #include <protocols/loops/loops_main.hh>

 #include <protocols/viewer/viewers.hh>

 #include <core/types.hh>

 #include <core/scoring/sasa.hh>

// #include <basic/prof.hh> // profiling
// #include <basic/CacheableData.hh> // profiling

 #include <core/id/SequenceMapping.hh>

 #include <core/chemical/AtomTypeSet.hh>
 #include <core/chemical/MMAtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>

 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/Ramachandran.hh>
 #include <core/pack/dunbrack/RotamerLibrary.hh>
 #include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
 #include <core/scoring/hbonds/HBondSet.hh>
 #include <core/scoring/hbonds/hbonds.hh>
 #include <core/scoring/etable/count_pair/CountPairFunction.hh>

 #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskOperation.hh>

 #include <core/kinematics/FoldTree.hh>
 #include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>
 #include <core/id/AtomID_Map.hh>

 #include <core/mm/MMTorsionLibrary.hh>
 #include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDB_PoseMap.hh>

#include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/PeriodicFunc.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <protocols/enzdes/EnzConstraintIO.hh>

// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>

//silly using/typedef


#include <basic/tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

 using basic::T;
 using basic::Error;
 using basic::Warning;

 static numeric::random::RandomGenerator RG(16621);

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
//using namespace protocols;

using utility::vector1;
using std::string;
using io::pdb::dump_pdb;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace pose;
using namespace chemical;
using namespace conformation;
using namespace scoring;
using namespace optimization;
using namespace kinematics;
using namespace protocols::moves;
using namespace id;
using namespace protocols::frags;

Size
get_n_pep_nbrs(
	pose::Pose const & pose,
	vector1< bool > const is_pep,
	Real const cutoff_cg
)
{
	Size n_pep_nbrs( 0 );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_pep[i] ) continue;
		bool cg_res_has_nbr( false );
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if( cg_res_has_nbr ) break;
			if ( is_pep[j] ) continue;
			Residue const & rsd2( pose.residue(j) );
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				if( cg_res_has_nbr ) break;
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) < cutoff_cg ) {
						cg_res_has_nbr = true;
						break;
					}
				}
			}
		}
		if( cg_res_has_nbr ) ++n_pep_nbrs;
	}
	return n_pep_nbrs;
}



void
dump_efactor_pdb(
	pose::Pose & pose,
	core::scoring::ScoreFunctionOP const scorefxn,
	std::string const & tag
)
{
	Size const nres( pose.total_residue() );
//	id::AtomID_Mask const & mask;
//	id::initialize( mask, pose );

	( *scorefxn )( pose );
	scorefxn->accumulate_residue_total_energies( pose );

	char const *filename = tag.c_str();
	std::fstream out( filename, std::ios::out );


	Size number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	out << "MODEL     " << tag << "\n";
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		Real residue_total_energy( pose.energies().residue_total_energies( i ).dot( scorefxn->weights() ) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

 //			if ( ! mask[ id::AtomID( j,i ) ] ) continue;

			//skip outputing virtual atom unless specified
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&
				(int)atom.type() == (int)rsd.atom_type_set().n_atomtypes() ) continue;

			++number;
			assert( rsd.chain() < int(chains.size()) ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2, residue_total_energy ) << '\n';
		}
	}
	out << "ENDMDL\n";
}

Real
GetRmsd(
	pose::Pose const & ref_pose,
	pose::Pose const & new_pose,
	vector1< bool > const score_seqpos,
	bool const ca_only
)
{
	Real rmsd( 0 );
	Size natoms( 0 );
	if( !ca_only ){
		for( Size seqpos = 1; seqpos <= new_pose.total_residue(); ++seqpos ){
			if( !score_seqpos[ seqpos ] ) continue;
			for( Size atomno = 1; atomno <= ref_pose.residue( seqpos ).natoms(); ++atomno ){
				rmsd += ref_pose.residue( seqpos ).xyz( atomno ).distance_squared(
						new_pose.residue( seqpos ).xyz( atomno ) );
				++natoms;
			}
		}
	}
	else{
		for( Size seqpos = 1; seqpos <= new_pose.total_residue(); ++seqpos ){
			if( !( score_seqpos[ seqpos ] && new_pose.residue( seqpos ).is_protein() ) ) continue;
			rmsd += ref_pose.residue( seqpos ).xyz( "CA" ).distance_squared(
					new_pose.residue( seqpos ).xyz( "CA" ) );
			++natoms;
		}
	}
	rmsd = std::sqrt( rmsd / natoms );
	return rmsd;
}


void
gen_fold_tree_for_nbr_segments(
	pose::Pose & pose,
	FoldTree & ftree,
	vector1< bool > const & is_ligand,
	vector1< bool > const & is_skipped,
	Real const & nbr_cutoff,
	vector1< bool > & is_nbr
)
{
	Size nres( pose.total_residue() );
	//define neighbors, all protein residues within cutoff from ligand, excluding is_skipped
	set_ss_from_phipsi( pose );
	for( Size i=1; i<= nres; ++i ) {
		if ( !is_ligand[ i ] ) continue;
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= nres; ++j ) {
			Residue const & rsd2( pose.residue(j) );
			if( is_ligand[ j ] || is_skipped[ j ] || !rsd2.is_protein() || !( pose.secstruct( j ) == 'L' ) ) continue;
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				if( is_nbr[ j ] ) break;
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < nbr_cutoff*nbr_cutoff ) {
						is_nbr[ j ] = true;
						break;
					}
				}
			}
		}
	}

	for( Size i = 2, j = 1; i <= nres; ++i, ++j ){
		if( is_nbr[ i ] && !is_nbr[ i-1 ] ) j = 1;
		if( is_nbr[ i ] && !is_nbr[ i+1 ] ){
			ftree.new_jump( i - j, i - static_cast< int >( j / 2 ), i - j );
			ftree.new_jump( i - j, i + 1, i );
		}
	}

}

bool
has_clash(
		pose::Pose & pose,
		vector1< bool > is_checked,
		scoring::ScoreFunctionOP const & scorefxn,
		Real const clash_threshold
	 )
{
	using namespace scoring;
	using namespace chemical;

	using namespace ObjexxFCL::format; // I and F

	bool is_clash( false );
	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( !is_checked[ seqpos ] ) continue;

		( *scorefxn )( pose );
		scorefxn->accumulate_residue_total_energies( pose );


		// cached energies object
		Energies & energies( pose.energies() );

		// the neighbor/energy links
		EnergyGraph & energy_graph( energies.energy_graph() );

		// search upstream
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node( seqpos )->lower_edge_list_begin(),
				irue = energy_graph.get_node( seqpos )->lower_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
			Size const j( edge->get_first_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
//				std::cout<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << std::endl;
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;

		// and downstream
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node( seqpos )->upper_edge_list_begin(),
				irue = energy_graph.get_node( seqpos )->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
//				std::cout<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << std::endl;
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;
	}
	return is_clash;
}




void
RunPepSpec()
{
	//load loop lengths//
	Size n_peptides( option[ pep_spec::n_peptides ] );
	Size n_build_loop( option[ pep_spec::n_build_loop ] );
	Size n_cgrelax_loop( option[ pep_spec::n_cgrelax_loop ]  );
	Size n_ramp_loop;
	if( option[ pep_spec::ramp ] ) n_ramp_loop = option[ pep_spec::n_ramp_loop ];
	Size n_farelax_loop( option[ pep_spec::n_farelax_loop ]  );

	//load input structure pdbs/csts//
	pose::Pose pose;
	std::string input_name;
	vector1< Pose > poses;
	vector1< std::string > pdb_filenames;
	vector1< std::string > cst_filenames;
	if( option[ pep_spec::pdb_list ].user() ){
		std::string pdb_list_filename( option[ pep_spec::pdb_list ] );
		std::ifstream pdb_list_data( pdb_list_filename.c_str() );
		if ( !pdb_list_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_filename + '\n' );
		}
		std::string pdb_list_line;
		while( !getline( pdb_list_data, pdb_list_line, '\n' ).eof() ) {
			std::string this_filename( pdb_list_line );
			pdb_filenames.push_back( this_filename );
			Pose this_pose;
			core::import_pose::pose_from_pdb( this_pose, this_filename );
			poses.push_back( this_pose );
		}
		input_name = pdb_filenames[ 1 ];
		if( option[ pep_spec::cst_list ].user() ){
			std::string cst_list_filename( option[ pep_spec::cst_list ] );
			std::ifstream cst_list_data( cst_list_filename.c_str() );
			if( !cst_list_data.good() ) {
				utility_exit_with_message( "Unable to open file: " + cst_list_filename + '\n' );
			}
			std::string cst_list_line;
			while( !getline( cst_list_data, cst_list_line, '\n' ).eof() ) {
				std::string this_cst_filename( cst_list_line );
				cst_filenames.push_back( this_cst_filename );
			}
			if( cst_filenames.size() != pdb_filenames.size() ) utility_exit_with_message( "Error in cst_list: Must have one cst filename per pdb filename\n" );
		}
	}
	else{
		input_name = basic::options::start_file();
		pdb_filenames.push_back( input_name );
		Pose this_pose;
		core::import_pose::pose_from_pdb( this_pose, input_name );
		poses.push_back( this_pose );

	}
	if( option[ OptionKeys::constraints::cst_file ].user() ){
		cst_filenames.push_back( option[ OptionKeys::constraints::cst_file ] );
	}

	core::import_pose::pose_from_pdb( pose, input_name );
	pose::Pose start_pose( pose );

	pose::Pose ref_pose;
	if( option[ pep_spec::rmsd ] ){
		std::string ref_name( option[ pep_spec::rmsd_ref ] );
		core::import_pose::pose_from_pdb( ref_pose, ref_name );
	}

	//convert user input to internal values//
        Size const nres( pose.total_residue() );
	Size const prot_anchor_in( option[ pep_spec::prot_anchor ] );
       	Size const prot_begin_in ( option[ pep_spec::prot_begin ] );
       	Size const prot_end_in ( option[ pep_spec::prot_end ] );
	std::string const prot_chain_in( option[ pep_spec::prot_chain ] );
	Size const pep_anchor_in( option[ pep_spec::pep_anchor ] );
       	Size const pep_begin_in ( option[ pep_spec::pep_begin ] );
       	Size const pep_end_in ( option[ pep_spec::pep_end ] );
	std::string const pep_chain_in( option[ pep_spec::pep_chain ] );

	Size const prot_anchor( pose.pdb_info()->pdb2pose( prot_chain_in[0], prot_anchor_in ) ); // cat-ASP
	Size const prot_begin( pose.pdb_info()->pdb2pose( prot_chain_in[0], prot_begin_in ) );
	Size const prot_end( pose.pdb_info()->pdb2pose( prot_chain_in[0], prot_end_in  ) );
	Size const pep_anchor( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_anchor_in ) );
	Size const pep_begin( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_begin_in ) );
	Size const pep_end( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_end_in ) );
     	Size const nres_pep( pep_end - pep_begin + 1 );

	Size ts_bond_seqpos_in;
	std::string ts_bond_chain_in;
	Size ts_bond_seqpos;
	if( option[ pep_spec::ts_bond ] ){
		ts_bond_seqpos_in = option[ pep_spec::ts_bond_seqpos ];
		ts_bond_chain_in = option[ pep_spec::ts_bond_chain ];
		ts_bond_seqpos = pose.pdb_info()->pdb2pose( ts_bond_chain_in[0], ts_bond_seqpos_in );
	}

	Real cutoff( option[ pep_spec::interface_cutoff ] );
	std::string output_seq;

//	dump_pose_kinemage( "qq.kin", pose );

	//define prot, pep, and anchor residues//
	vector1< bool > is_pep( nres, false ), is_pep_not_anchor( nres, false ), is_prot( nres, false ), is_anchor( nres, false );
	for ( Size i=1; i<= nres; ++i ) {
		is_pep[i] = ( i >= pep_begin && i <= pep_end );
		is_prot[i] = ( i >= prot_begin && i <= prot_end );
		is_anchor[i] = ( (i==prot_anchor) || (i==pep_anchor) );
		is_pep_not_anchor[i] = ( i >= pep_begin && i <= pep_end && !(i==pep_anchor) );
	}

	//gen fold tree//
	FoldTree f( pose.total_residue() );

	vector1< bool > is_flex_prot( nres, false );
	if( option[ pep_spec::flex_prot_bb ] ){
		gen_fold_tree_for_nbr_segments( pose, f, is_pep, is_anchor, cutoff, is_flex_prot );
		std::cout << "sel flex_prot, ";
		for( Size i = 1; i <= nres; ++i ){
			if( is_flex_prot[ i ] ) std::cout << "resi " << string_of( i ) << " or ";
		}
		std::cout << std::endl;
	}

	if( pose.chain( prot_anchor ) < pose.chain( pep_anchor ) || ( pose.chain( prot_anchor ) == pose.chain( pep_anchor ) && prot_end < pep_begin ) ){
		Size pep_jump( f.new_jump( prot_anchor, pep_anchor, pep_begin - 1 ) );
		f.set_jump_atoms( pep_jump, "CA", "CB" );
		if( pep_end != pose.total_residue() ) f.new_jump( prot_anchor, pep_end + 1, pep_end );
		if( prot_end + 1 != pep_begin ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
	}
	else{
		Size pep_jump( f.new_jump( pep_anchor, prot_anchor, prot_begin - 1 ) );
		f.set_jump_atoms( pep_jump, "CB", "CA" );
		if( prot_end != pose.total_residue() ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
		if( pep_end + 1 != prot_begin ) f.new_jump( pep_anchor, pep_end + 1, pep_end );
	}

	//define scoring functions//
	core::scoring::ScoreFunctionOP soft_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::soft_wts ] ) );
	core::scoring::ScoreFunctionOP full_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
	core::scoring::ScoreFunctionOP cen_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::cen_wts ] ) );
	core::scoring::ScoreFunctionOP final_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::final_wts ] ) );

	Real original_full_fa_rep( full_scorefxn->get_weight( fa_rep ) );

	//can probably delete this ramp stuff//
	Real original_full_fa_atr( full_scorefxn->get_weight( fa_atr ) );
	Real original_soft_fa_rep = soft_scorefxn->get_weight( fa_rep );
	Real original_soft_fa_atr = soft_scorefxn->get_weight( fa_atr );
	Real min_full_fa_rep;
	Real min_full_fa_atr;
	Real min_soft_fa_rep;
	Real min_soft_fa_atr;
	if( option[ pep_spec::ramp ] ){
		min_full_fa_rep = original_full_fa_rep * option[ pep_spec::min_fa_rep_frac ];
		min_full_fa_atr = original_full_fa_atr * option[ pep_spec::min_fa_rep_frac ];
		min_soft_fa_rep = original_soft_fa_rep * option[ pep_spec::min_fa_rep_frac ];
		min_soft_fa_atr = original_soft_fa_atr * option[ pep_spec::min_fa_rep_frac ];
	}

	//score prot ref state//
	Pose prot_eval_pose( pose );
	Real prot_score;
	if( option[ pep_spec::score_binding ] ){
		Size eval_prot_begin( prot_begin );
		Size eval_prot_end( prot_end );
		if( pep_end < prot_begin ){
			eval_prot_begin -= ( pep_end - pep_begin + 1 );
			eval_prot_end -= ( pep_end - pep_begin + 1 );
		}
		prot_eval_pose.conformation().delete_residue_range_slow( pep_begin, pep_end );

		(*final_scorefxn)( prot_eval_pose );
		final_scorefxn->accumulate_residue_total_energies( prot_eval_pose );
		prot_score = prot_eval_pose.energies().total_energies().dot( final_scorefxn->weights() );
/*
		if( option[ pep_spec::score_binding_repack ] ){
			pack::task::PackerTaskOP prot_eval_task( pack::task::TaskFactory::create_packer_task( prot_eval_pose ));
			prot_eval_task->initialize_from_command_line();
			for( Size seqpos = 1; seqpos <= prot_eval_pose.total_residue(); ++seqpos ){
				if( seqpos >= eval_prot_begin && seqpos <= eval_prot_end ) prot_eval_task->nonconst_residue_task( seqpos ).restrict_to_repacking();
				else prot_eval_task->nonconst_residue_task( seqpos ).prevent_repacking();
			}
			protocols::simple_moves::PackRotamersMoverOP prot_eval_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, prot_eval_task, 1 ) );
			prot_eval_pack->apply( prot_eval_pose );
		}
		if( option[ pep_spec::score_binding_min ] ){
			kinematics::MoveMapOP mm_prot_eval ( new kinematics::MoveMap );
			mm_prot_eval->set_chi( true );
			protocols::simple_moves::MinMoverOP prot_eval_min_mover = new protocols::simple_moves::MinMover( mm_prot_eval, full_scorefxn, "dfpmin", 0.001, true );
			prot_eval_min_mover->apply( prot_eval_pose );
			prot_score = prot_eval_pose.energies().total_energies().dot( full_scorefxn->weights() );
			dump_pdb( prot_eval_pose,  option[ out::file::o ] + "_prot.pdb" );
		}
		std::cout<<"!!prot "+ option[ out::file::o ] + "_prot.pdb" + "\t"<<prot_eval_pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() );
		std::cout<<"\ttotal_score:\t"<<prot_score<<std::endl;
*/
	}

	//////////////build cen pep/////////////////
	VallData vall( option[ pep_spec::vall ] );
	TorsionFragmentLibrary lib;
	Size const frag_size( 1 );

	lib.resize( nres_pep - frag_size + 1 );
	for ( Size i=1; i<= nres_pep - frag_size + 1; ++i ) {
		Size seqpos( i + pep_begin - 1 );

		std::string frag_seq;
		std::string ss_seq;
		for( Size ii = 1; ii <= frag_size; ++ii ){
			frag_seq.append( 1, pose.residue( seqpos ).name1() );
			ss_seq.append( 1, 'S' );
		}

		Real seq_wt( 0.0 );
//		if( option[ pep_spec::use_input_seq ] || seqpos == pep_anchor ){
		if( option[ pep_spec::use_input_seq ] ){
			seq_wt = 10.0;
		}
		vall.get_frags( 10000, frag_seq, ss_seq, seq_wt, 1.0, false, false, true, lib[i] );

	}


	for( Size peptide_loop = 1; peptide_loop <= n_peptides; ++peptide_loop ){	//LOOP how many final peps

		//load random start pdb and cst//
		int pose_index( static_cast< int >( RG.uniform() * poses.size() + 1 ) );
		std::cout<<"!! Initializing "<< option[ out::file::o ] + "_" + string_of( peptide_loop ) + " with " + pdb_filenames[ pose_index ] << std::endl;
		pose = poses[ pose_index ];
		std::string cst_filename;
		if( option[ pep_spec::cst_list ].user() ) cst_filename = cst_filenames[ pose_index ];
		if( option[ OptionKeys::constraints::cst_file ].user() ) cst_filename = cst_filenames[ 1 ];

		//set foldtree in the pose//
		pose.fold_tree( f );

		//use input sequence//
		if( option[ pep_spec::use_input_seq ] ){
			std::string const input_seq( option[ pep_spec::input_seq ] );
			for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
				if(mut_site==pep_anchor) continue;
				chemical::AA input_aa( aa_from_oneletter_code( input_seq[ mut_site - pep_begin ] ) );
				chemical::make_sequence_change( mut_site, input_aa, pose );
			}
		}
		//or erase identity//
		else if( !( option[ pep_spec::test_no_design ] ) ){
			for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
				if(mut_site==pep_anchor) continue;
				chemical::make_sequence_change( mut_site, chemical::aa_ala, pose );
			}
		}

		//convert to CG residues//
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		( *cen_scorefxn )( pose );
		cen_scorefxn->accumulate_residue_total_energies( pose );

		//make backbone from random phi/psi insertions//

		MonteCarloOP mc_frag ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );
		if( !option[ pep_spec::test_no_frag ] ){
			MonteCarloOP mc_frag ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );
			for ( Size build_loop = 1; build_loop <= static_cast< int >( std::sqrt( n_build_loop ) ); ++build_loop ) {
				mc_frag->reset( pose );
				for ( Size build_loop_inner = 1; build_loop_inner <= static_cast< int >( std::sqrt( n_build_loop ) ); ++build_loop_inner ) {
					// choose an insertion position
					Size pos(  static_cast< int > ( nres_pep * RG.uniform() ) + 1 );

					Size const nfrags( lib[ pos ].size() );
					int const frag_index( static_cast< int >( nfrags * RG.uniform() + 1 ) );
					lib[ pos ][ frag_index ].insert( pose, pep_begin - 1 + pos );

					mc_frag->boltzmann( pose );
				}

			}
			mc_frag->recover_low( pose );
		}

		////////////////????????????????????
		/*
		// Try to extend chain
		ResidueTypeSet const & cg_rsd_set( pose.residue(1).residue_type_set() );
		ResidueOP ala( ResidueFactory::create_residue( cg_rsd_set.name_map( "ALA" ) ) );
		pose.conformation().safely_append_polymer_residue_after_seqpos( *ala, pep_end, true );
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *ala, pep_begin, true );
		Size pep_end_cg( pep_end + 2 );
		Size pep_anchor_cg( pep_anchor + 1 );
		Size nres_pep_cg( nres_pep + 2 );
		pose.set_omega( pep_begin, 180.0 );
		pose.set_omega( pep_end_cg - 1, 180.0 );

		MonteCarloOP mc_extend ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );
		for ( Size build_loop = 1; build_loop <= n_build_loop; ++build_loop ) {
				// choose an insertion position
				Size pos(  static_cast< int > ( RG.uniform() + 1 ) * ( nres_pep_cg - 1 ) + 1 );

				Size const nfrags( lib[ 1 ].size() );
				int const frag_index( static_cast< int >( nfrags * RG.uniform() + 1 ) );
				lib[ 1 ][ frag_index ].insert( pose, pep_begin - 1 + pos );

				mc_extend->boltzmann( pose );
		}
		mc_extend->recover_low( pose );

		if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, cen_scorefxn, option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_precen.pdb" );
		std::cout<<"!!! "+ option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_precen.pdb" +
				"\t"<<pose.energies().total_energies().weighted_string_of( cen_scorefxn->weights() );
		std::cout<<"\ttotal_score:\t"<<pose.energies().total_energies()[ total_score ]<<"\n";

		pose.conformation().delete_residue_slow( pep_end_cg );
		pose.conformation().delete_residue_slow( pep_begin );
*/
////////////////????????????????????

		//debug dump CG data and pdb//
		if( option[ pep_spec::test_frag_only ] ){
			if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, cen_scorefxn, option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_cen.pdb" );
			Real cg_rmsd;
			if( option[ pep_spec::rmsd ] ){
				cg_rmsd = GetRmsd( ref_pose, pose, is_pep, true );
			}
			std::cout<<"!!! "+ option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_cen.pdb" +
				"\t"<<pose.energies().total_energies().weighted_string_of( cen_scorefxn->weights() );
			std::cout<<"\ttotal_score:\t"<<pose.energies().total_energies()[ total_score ];

			if( option[ pep_spec::rmsd ] ){
				std::cout<<"\tca_rmsd:\t"<< string_of( cg_rmsd );
			}
			std::cout<<std::endl;
			continue;
		}
		///////////////////////////////////////////////////////////////

		//make peptide from CG backbone//
		Pose restart_cgrelax_pose( pose );
		for(Size cgrelax_loop = 1; cgrelax_loop <= n_cgrelax_loop; cgrelax_loop++){
			pose = restart_cgrelax_pose;

			//min movemap//
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			if( !option[ pep_spec::test_no_min_bb ] ){
				mm_min->set_bb( is_pep );
			}
			if( option[ pep_spec::flex_prot_bb ] ) mm_min->set_chi( is_flex_prot );
			if( option[ pep_spec::flex_prot_bb ] ) mm_min->set_bb( is_flex_prot );

			//small-shear movemap//
			kinematics::MoveMapOP mm_move ( new kinematics::MoveMap );
			mm_move->set_chi( is_pep );
			if( !option[ pep_spec::test_no_min_bb ] ){
				mm_move->set_bb( is_pep );
				if( option[ pep_spec::flex_prot_bb ] ) mm_move->set_bb( is_flex_prot );
			}

			//small CG backbone moves//
			if( !( option[ pep_spec::test_no_cgrelax ] || cgrelax_loop == 1 ) ){

				protocols::simple_moves::SmallMoverOP cg_small( new protocols::simple_moves::SmallMover( mm_move, 2.0, 1 ) ); //LOOP
				cg_small->angle_max( 'H', 2.0 );
				cg_small->angle_max( 'E', 2.0 );
				cg_small->angle_max( 'L', 2.0 );

				MonteCarloOP mc_cg ( new MonteCarlo( pose, *cen_scorefxn, 2.0  ) );
				TrialMoverOP cg_trial = new TrialMover( cg_small, mc_cg );
				RepeatMoverOP cg_cycle = new RepeatMover( cg_trial, 100 );

				for( Size cgmove_loop = 1; cgmove_loop <= 10; cgmove_loop++ ){
					mc_cg->reset( pose );
					cg_cycle->apply( pose );
				}
				mc_cg->recover_low( pose );
			}
			//switch back to fullatom
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );

			//replace kin residues w/ original rotamers
			for(Size resnum = prot_begin; resnum <= prot_end; ++resnum){
				pose.replace_residue( resnum, start_pose.residue( resnum ), false );
			}
			pose.replace_residue( pep_anchor, start_pose.residue( pep_anchor ), true );

			if( option[ pep_spec::ts_bond ] ){
				core::pose::add_variant_type_to_pose_residue( pose, "ATP_PG_CONNECT", ts_bond_seqpos );
				core::pose::add_variant_type_to_pose_residue( pose, "SER_OG_CONNECT", pep_anchor );
				pose.conformation().declare_chemical_bond( ts_bond_seqpos, "PG", pep_anchor, "OG" );
			}

			//turn on constraints//
			using namespace core::scoring::constraints;
			if( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::pep_spec::cst_list ].user() ){
				soft_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );
				full_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );
				ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( cst_filename, new ConstraintSet, pose );
				pose.constraint_set( cst_set );
			}

			//define neighbors
			vector1< bool > is_pep_nbr( nres, false );
			for ( Size i=1; i<= nres; ++i ) {
				Residue const & rsd1( pose.residue(i) );
				if ( is_pep[i] ) continue;
				for ( Size j=1; j<= nres; ++j ) {
					Residue const & rsd2( pose.residue(j) );
					if ( !is_pep[j] ) continue;
					if( is_pep_nbr[i] ) break;
					for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
						if( is_pep_nbr[i] ) break;
						for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
							if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) < cutoff ) {
								is_pep_nbr[i] = true;
								break;
							}
						}
					}
				}
			}
			mm_min->set_chi( is_pep_nbr );

			//define design task and repack task
			vector1< bool > design_this( nres, false );
			vector1< bool > repack_this( nres, false );
			vector1< bool > check_clash( nres, true );
			pack::task::PackerTaskOP dz_task( pack::task::TaskFactory::create_packer_task( pose ));
			dz_task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= nres; ++i ) {
				if ( is_pep[ i ] && !( is_anchor[i] ) ) {
					if( ( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ){
						repack_this[ i ] = true;
						dz_task->nonconst_residue_task( i ).restrict_to_repacking();
					}
					else{
						design_this[ i ] = true;
					}
				} else if ( is_pep_nbr[i] && is_prot[i] && !( is_anchor[i]) ) {
					repack_this[ i ] = true;
					dz_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					dz_task->nonconst_residue_task( i ).prevent_repacking();
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( new pack::task::DesignOrRepackOperation( design_this, repack_this ) );

/*

			//poly-A FA refinement//
			if( !( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ){
				for( Size mut_site = 1; mut_site <= nres; mut_site++ ){
					if( is_anchor[ mut_site ] || !( is_pep[ mut_site ] || is_prot[ mut_site ] ) ) continue;
					if( is_pep[ mut_site ] ){
						chemical::make_sequence_change( mut_site, chemical::aa_gly, pose );
					}
					else if( is_pep_nbr[ mut_site ] && pose.residue( mut_site ).name1() != 'G' ){
						chemical::make_sequence_change( mut_site, chemical::aa_ala, pose );
					}
				}
			}
			core::scoring::ScoreFunctionOP rep_scorefxn( new ScoreFunction );
			rep_scorefxn->set_weight( fa_rep, original_full_fa_rep );
			rep_scorefxn->set_weight( rama, 1.0 );
			( *rep_scorefxn )( pose );
			rep_scorefxn->accumulate_residue_total_energies( pose );

			if( has_clash( pose, is_pep_not_anchor, rep_scorefxn, 5.0 ) ){
				protocols::simple_moves::SmallMoverOP rep_small_mover( new protocols::simple_moves::SmallMover( mm_move, 1.0, 1 ) );	//LOOP
				rep_small_mover->angle_max( 'H', 5.0 );
				rep_small_mover->angle_max( 'E', 5.0 );
				rep_small_mover->angle_max( 'L', 5.0 );
				MonteCarloOP mc_rep ( new MonteCarlo( pose, *rep_scorefxn, 1.0 ) );
				TrialMoverOP rep_trial = new TrialMover( rep_small_mover, mc_rep );
				RepeatMoverOP rep_cycle = new RepeatMover( rep_trial, 1000 );

				rep_cycle->apply( pose );
				mc_rep->recover_low( pose );
			}
			if( has_clash( pose, is_pep_not_anchor, rep_scorefxn, 5.0 ) ){
				protocols::simple_moves::MinMoverOP rep_min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );
				rep_min_mover->apply( pose );
			}
*/

			//replace kin residues w/ original rotamers
			for(Size resnum = prot_begin; resnum <= prot_end; ++resnum){
				pose.replace_residue( resnum, start_pose.residue( resnum ), false );
			}

/*
			if( option[ pep_spec::ramp ] && RG.uniform() < option[ pep_spec::ramp_prob ] ){
				std::cout << "!! Ramping fa_rep for " << option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + ".pdb" << std::endl;
				if( option[ pep_spec::test_dump_all_ramp ] ) dump_efactor_pdb( pose, full_scorefxn,  option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_preramp.pdb" );
//				core::scoring::ScoreFunctionOP clash_scorefxn( full_scorefxn );
//				core::scoring::ScoreFunctionOP soft_clash_scorefxn( soft_scorefxn );
				core::scoring::ScoreFunctionOP clash_scorefxn( new ScoreFunction );
				core::scoring::ScoreFunctionOP soft_clash_scorefxn( new ScoreFunction );
				soft_clash_scorefxn->set_etable( FA_STANDARD_SOFT );
				for( Size ramp_loop = 1; ramp_loop <= n_ramp_loop + 1; ++ramp_loop ){

					protocols::simple_moves::PackRotamersMoverOP ramp_pack( new protocols::simple_moves::PackRotamersMover( soft_clash_scorefxn, dz_task, 1 ) );
					( *soft_clash_scorefxn )( pose );
					soft_clash_scorefxn->accumulate_residue_total_energies( pose );
					if( option[ pep_spec::use_input_seq ] )ramp_pack->apply( pose );

					Real ramp_factor( pow( static_cast< Real >( ramp_loop - 1 ) / static_cast< Real >( n_ramp_loop ), option[ pep_spec::ramp_exp ] ) );
					Real new_full_fa_rep ( min_full_fa_rep + ( original_full_fa_rep - min_full_fa_rep ) * ramp_factor );
					Real new_full_fa_atr ( min_full_fa_atr + ( original_full_fa_atr - min_full_fa_atr ) * ramp_factor );
					Real new_soft_fa_rep ( min_soft_fa_rep + ( original_soft_fa_rep - min_soft_fa_rep ) * ramp_factor );
					Real new_soft_fa_atr ( min_soft_fa_atr + ( original_soft_fa_atr - min_soft_fa_atr ) * ramp_factor );
					clash_scorefxn->set_weight( fa_rep, new_full_fa_rep  );
					clash_scorefxn->set_weight( fa_atr, new_full_fa_atr  );
					soft_clash_scorefxn->set_weight( fa_rep, new_soft_fa_rep  );
					soft_clash_scorefxn->set_weight( fa_atr, new_soft_fa_atr  );
					protocols::simple_moves::MinMoverOP clash_min_mover = new protocols::simple_moves::MinMover( mm_min, soft_clash_scorefxn, "dfpmin", 0.001, true );
					clash_min_mover->apply( pose );
					if( option[ pep_spec::test_dump_all_ramp ] ) dump_efactor_pdb( pose, full_scorefxn,  option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_ramp" + string_of( ramp_loop - 1 ) + ".pdb" );
				}
			}
			else{
				std::cout << "!! Not ramping fa_rep for " << option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) << std::endl;
			}
*/

			//randomize pep sequence//
			if( !( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ){
				for(Size mut_site = pep_begin; mut_site <= pep_end; mut_site++){ //over all pep positions
					if(mut_site==pep_anchor) continue;
					int resindex;
					resindex = static_cast< int > ( 20 * RG.uniform() + 1 );
					chemical::make_sequence_change( mut_site, chemical::AA(resindex), pose );
				}
			}

			//allow jump relaxation//
			if( !option[ pep_spec::freeze_pep_anchor ] ) mm_min->set_jump( 1, true );

			//define movers//
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

			protocols::simple_moves::ShearMoverOP shear_mover( new protocols::simple_moves::ShearMover( mm_move, 5.0, 1 ) );	//LOOP
			shear_mover->angle_max( 'H', 5.0 );
			shear_mover->angle_max( 'E', 5.0 );
			shear_mover->angle_max( 'L', 5.0 );

			protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, dz_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				design_seq->add_mover( dz_pack );
				design_seq->add_mover( dz_rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				design_seq->add_mover( min_mover );
			}

			SequenceMoverOP relax_seq = new SequenceMover;
			relax_seq->add_mover( shear_mover );
			relax_seq->add_mover( dz_rottrial );
			if( !option[ pep_spec::test_no_min ] ){
				relax_seq->add_mover( min_mover );
			}

			//RUN MOVERS//
			( *soft_scorefxn )( pose );
			soft_scorefxn->accumulate_residue_total_energies( pose );
			design_seq->apply( pose );

			MonteCarloOP mc_relax ( new MonteCarlo( pose, *full_scorefxn, 1.0 ) );
			TrialMoverOP relax_trial = new TrialMover( relax_seq, mc_relax );
			RepeatMoverOP relax_cycle = new RepeatMover( relax_trial, n_farelax_loop );
			if( !option[ pep_spec::test_no_farelax ] ){
				relax_cycle->apply( pose );
			}
			mc_relax->recover_low( pose );

			( *final_scorefxn )( pose );
			final_scorefxn->accumulate_residue_total_energies( pose );

			//Analysis//
			output_seq.clear();
			for(Size i = pep_begin; i <= pep_end; i++){
				output_seq.append( 1, pose.residue( i ).name1() );
			}
			Real total_score( pose.energies().total_energies().dot( final_scorefxn->weights() ) );

			Real ca_rmsd;
			Real fa_rmsd;
			if( option[ pep_spec::rmsd ] ){
				ca_rmsd = GetRmsd( ref_pose, pose, is_pep, true );
				if( ( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ) fa_rmsd = GetRmsd( ref_pose, pose, is_pep, false );
			}


			//score pep unbound state//
			Pose pep_eval_pose( pose );
			Real pep_score;
			Real bind_score;
			bool abort( false );
			if( option[ pep_spec::score_binding ] ){
				//remove ts bond//
				if( option[ pep_spec::ts_bond ] ){
					core::pose::remove_variant_type_from_pose_residue( pep_eval_pose, "SER_OG_CONNECT", pep_anchor );
				}
				//remove prot//
				if( pep_begin != 1 ){
					pep_eval_pose.conformation().delete_residue_range_slow( 1, pep_begin - 1 );
					if( pep_eval_pose.conformation().chain_end( 1 ) != pep_eval_pose.total_residue() ){
						pep_eval_pose.conformation().delete_residue_range_slow( pep_eval_pose.conformation().chain_end( 1 ) + 1, pep_eval_pose.total_residue() );
					}
				}
				else if( pep_end != pep_eval_pose.total_residue() ){
					pep_eval_pose.conformation().delete_residue_range_slow( pep_end + 1, pep_eval_pose.total_residue() );
				}
				//calc unbound score//
				Real best_pep_score( 0.0 );
				Pose best_pep_eval_pose( pep_eval_pose );
				Pose pep_restart_pose = pep_eval_pose;
				Size n_unbound_loop( option[ pep_spec::n_unbound_loop ] );
				for( Size unbound_loop = 1; unbound_loop <= n_unbound_loop; ++unbound_loop ){
					pep_eval_pose = pep_restart_pose;
					if( option[ pep_spec::score_binding_refold ] && unbound_loop != 1 ){
						core::util::switch_to_residue_type_set( pep_eval_pose, core::chemical::CENTROID );
						core::scoring::ScoreFunctionOP pep_cen_scorefxn(  ScoreFunctionFactory::create_score_function( "cen_std.wts" ) );
						( *pep_cen_scorefxn )( pep_eval_pose );
						pep_cen_scorefxn->accumulate_residue_total_energies( pep_eval_pose );

						vector1< core::scoring::Ramachandran > pep_rama_movers;
						core::scoring::Ramachandran loop_pep_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
						loop_pep_rama_mover.init_rama_sampling_table_for_ss_type( 3 );
						pep_rama_movers.push_back( loop_pep_rama_mover );
						if( option[ pep_spec::random_rama_ss_type ] ){
							core::scoring::Ramachandran helix_pep_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
							helix_pep_rama_mover.init_rama_sampling_table_for_ss_type( 1 );
							pep_rama_movers.push_back( helix_pep_rama_mover );
							core::scoring::Ramachandran sheet_pep_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
							sheet_pep_rama_mover.init_rama_sampling_table_for_ss_type( 2 );
							pep_rama_movers.push_back( sheet_pep_rama_mover );
						}
						MonteCarloOP mc_pep_cen ( new MonteCarlo( pep_eval_pose, *pep_cen_scorefxn, 2.0 ) );
						mc_pep_cen->reset( pep_eval_pose );
						for ( Size pep_build_loop = 1; pep_build_loop <= n_build_loop; ++pep_build_loop ){
							Size pep_seqpos(  static_cast< int > ( nres_pep * RG.uniform() ) + 1 );
							Size pep_rama_mover_index( 1 );
							if( option[ pep_spec::random_rama_ss_type ] ) pep_rama_mover_index = static_cast< int > ( RG.uniform() * pep_rama_movers.size() + 1 );
							Real pep_rama_phi, pep_rama_psi;
							chemical::AA pep_aa( aa_from_oneletter_code( pep_eval_pose.residue( pep_seqpos ).name1() ) );
							pep_rama_movers[ pep_rama_mover_index ].random_phipsi_from_rama( pep_aa, pep_rama_phi, pep_rama_psi );
							pep_eval_pose.set_phi( pep_seqpos, pep_rama_phi );
							pep_eval_pose.set_psi( pep_seqpos, pep_rama_psi );
							mc_pep_cen->boltzmann( pep_eval_pose );
						}
						mc_pep_cen->recover_low( pep_eval_pose );
						core::util::switch_to_residue_type_set( pep_eval_pose, core::chemical::FA_STANDARD );
					}
					if( option[ pep_spec::score_binding_repack ] ){
						pack::task::PackerTaskOP pep_eval_task( pack::task::TaskFactory::create_packer_task( pep_eval_pose ));
						pep_eval_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
						protocols::simple_moves::PackRotamersMoverOP pep_eval_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, pep_eval_task, 1 ) );
						pep_eval_pack->apply( pep_eval_pose );
					}
					if( option[ pep_spec::score_binding_min ] ){
						( *final_scorefxn )( pep_eval_pose );
						final_scorefxn->accumulate_residue_total_energies( pep_eval_pose );
						kinematics::MoveMapOP mm_pep_eval ( new kinematics::MoveMap );
						mm_pep_eval->set_chi( true );
//						mm_pep_eval->set_bb( true );
						protocols::simple_moves::MinMoverOP pep_eval_min_mover = new protocols::simple_moves::MinMover( mm_pep_eval, final_scorefxn, "dfpmin", 0.001, true );
						Real premin_pep_score = pep_eval_pose.energies().total_energies().dot( final_scorefxn->weights() );
						//Pose pep_eval_premin_pose( pep_eval_pose );
						pep_eval_min_mover->apply( pep_eval_pose );
						pep_score = pep_eval_pose.energies().total_energies().dot( final_scorefxn->weights() );
						// ignore if minimizer fails
						if( pep_score > premin_pep_score ) abort = true;
					}
					( *final_scorefxn )( pep_eval_pose );
					final_scorefxn->accumulate_residue_total_energies( pep_eval_pose );
					pep_score = pep_eval_pose.energies().total_energies().dot( final_scorefxn->weights() );
					if( pep_score < best_pep_score || unbound_loop == 1 ){
						best_pep_score = pep_score;
						best_pep_eval_pose = pep_eval_pose;
					}
				}
				pep_score = best_pep_score;
				pep_eval_pose = best_pep_eval_pose;
//				bind_score = total_score - prot_score - pep_score;
				bind_score = total_score - pep_score;
			}
			if( abort ){
				std::cout << "!! min_fail\n";
				continue;
			}

			//Output//
			std::cout<<"!!! "+ option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + ".pdb" + "\t" << output_seq << "\t" << pose.energies().total_energies().weighted_string_of( final_scorefxn->weights() );
			std::cout<<"\ttotal_score:\t"<<total_score;
			if( option[ pep_spec::score_binding ] ){
				std::cout<<"\tbinding_score:\t"<<bind_score;
			}
			if( option[ pep_spec::rmsd ] ){
				std::cout<<"\tca_rmsd:\t"<< string_of( ca_rmsd );
				if( ( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ) std::cout<<"\tfa_rmsd:\t"<< string_of( fa_rmsd );
			}
			std::cout<<std::endl;

			if( option[ pep_spec::score_binding ] ){
				std::cout<<"!!pep "+ option[ out::file::o ] + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + "_pep.pdb" + "\t"<<output_seq<<"\t"<<pep_eval_pose.energies().total_energies().weighted_string_of( final_scorefxn->weights() );
				std::cout<<"\ttotal_score:\t"<<pep_score<<std::endl;
			}

			if( option[ pep_spec::test_dump_all ] ){
				dump_efactor_pdb( pose, full_scorefxn, option[ out::file::o ] + "_" + string_of( peptide_loop )
						+ "_" + string_of( cgrelax_loop ) + ".pdb" );
			}

			std::cout<<"\n";
		}
	}
}


int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		RunPepSpec();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
