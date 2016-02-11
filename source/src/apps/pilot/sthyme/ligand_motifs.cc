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
//#include <basic/options/option.hh>

#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/residue_io.hh>

#include <core/chemical/VariantType.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>


#include <core/graph/Graph.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/dna/util.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/Motif.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace ObjexxFCL;
using namespace pose;
using namespace chemical;
using namespace scoring;
using namespace optimization;

using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
output_single_motif(
	Pose & src_pose,
	std::string & pdb_name,
	int prot_pos,
	utility::vector1< Size > &  contacts
)
{

	int prot_pos_pdb( src_pose.pdb_info()->number( prot_pos ) );
	char prot_pos_chain( src_pose.pdb_info()->chain( prot_pos ) );

	std::string output_path( "Ligand_motif_dir/" );
	std::string delimiter( "_" );
	std::string extension( ".pdb" );
	std::string motif_file_name(
		output_path +
		src_pose.residue_type( prot_pos ).name3() + string_of( prot_pos_pdb ) + string_of( prot_pos_chain ) + delimiter );

	for( Size ic(1) ; ic <= contacts.size() ; ++ic ) {
		protocols::motifs::Motif motif( src_pose.residue( prot_pos ), src_pose.residue( contacts[ic] ) );
		conformation::ResidueOP refdnares = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
		conformation::Residue protres( src_pose.residue( prot_pos ) );
		conformation::Residue dnares( src_pose.residue( contacts[ic] ) );
		motif.place_residue( *refdnares, protres );
		motif.place_residue( protres, dnares );
		Pose pose;
		pose.append_residue_by_jump( protres, 1 );
		//pose.append_residue_by_jump( src_pose.residue( prot_pos ), 1 );
		Size dna_pos = contacts[ ic ];
		pose.append_residue_by_jump( dnares, 1 );
		pose.conformation().insert_chain_ending( 1 );
		int dna_pos_pdb( src_pose.pdb_info()->number( dna_pos ) );
		char dna_pos_chain( src_pose.pdb_info()->chain( dna_pos ) );
		std::string motif_file_name_2 = motif_file_name + src_pose.residue_type( dna_pos ).name1() + string_of( dna_pos_pdb ) + string_of( dna_pos_chain ) + delimiter;

		std::string motif_file_name_3 = motif_file_name_2 + pdb_name + extension;

		std::cout << "Writing " << motif_file_name_3 << std::endl;
		io::pdb::dump_pdb( pose, motif_file_name_3 );
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
place_waters_and_minimize( Pose & pose )
{

  ScoreFunction scorefxn;
  // Jim had these weights in his original motif collecting app - I'll leave them
	scorefxn.set_weight( fa_atr, 1.0806 );
  scorefxn.set_weight( fa_rep, 0.65 );
  scorefxn.set_weight( fa_sol, 0.5842 );
  scorefxn.set_weight( fa_dun, 0.7506 );

  scorefxn.set_weight( h2o_intra, 0.01 );
  scorefxn.set_weight( h2o_hbond, 1.0 );

  scorefxn.set_weight( hbond_sr_bb, 3.1097 );
  scorefxn.set_weight( hbond_lr_bb, 3.1097 );
  scorefxn.set_weight( hbond_bb_sc, 3.1097 );
  scorefxn.set_weight( hbond_sc, 3.1097 );

	Energy score_orig = scorefxn( pose );

	std::cout << "Score before pack " << score_orig << std::endl;

	// The packing is strictly to allow waters to be re-placed into the model
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	kinematics::MoveMap mm;
  mm.set_bb( false );
  mm.set_chi( false );
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if ( pose.residue(i).is_protein() ) {
			task->nonconst_residue_task(i).prevent_repacking();
			// only allowing protein minimization because I saw strange DNA changes when I allowed everything to minimize
			//mm.set_chi( i, true );
		}
	}
	task->set_bump_check( true );
	pack::pack_rotamers( pose, scorefxn, task);

	Energy end_score = scorefxn( pose );
	std::cout << "Score after pack " << end_score << std::endl;
	//io::pdb::dump_pdb( pose, "post_pack.pdb" );

	// Minimization depends on the weight set used, so therefore I am nervous to use it at all, in case it disturbs native interactions
	//MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1.0e-3, true /*use_nblist*/, false /*deriv_check*/ );

	//AtomTreeMinimizer minimizer;
	//minimizer.run( pose, mm, scorefxn, options );
  //io::pdb::dump_pdb( pose, "post_minimization.pdb" );

	//Energy min_score = scorefxn( pose );
	//std::cout << "Score after minimization " << min_score << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Real
get_packing_score(
	Pose & pose,
	Size pos1,
	Size pos2
)
{

	ScoreFunction pack_score;
  pack_score.set_weight( fa_atr, 0.80 );
  pack_score.set_weight( fa_rep, 0.44 );

	//TwoBodyEnergyMap pack_map;

	//pack_score.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, pack_map );

	//return pack_map[ fa_atr ];
}

Real
get_hbond_score(
	Pose & pose,
	Size pos1,
	Size pos2
)
{

	ScoreFunction hbond_score;
  hbond_score.set_weight( hbond_sc, 1.0 );

	//TwoBodyEnergyMap hbond_map;

	//hbond_score.eval_cd_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, hbond_map );

	//return hbond_map[ hbond_sc ];
}

Real
get_water_hbond_score(
	Pose & pose,
	Size pos1,
	Size pos2
)
{

	ScoreFunction water_hbond_score;
  water_hbond_score.set_weight( h2o_hbond, 1.0 );

	//TwoBodyEnergyMap water_hbond_map;

	//water_hbond_score.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, water_hbond_map );

	//return water_hbond_map[ h2o_hbond ];
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
process_for_motifs(
	Pose & pose,
	std::string & pdb_name,
	protocols::motifs::MotifLibrary & motifs
)
{
	int nres( pose.total_residue() );

	// Place waters and minimize
	//place_waters_and_minimize( pose );

	// Loop over positions, skipping non-amino acid
	for( int prot_pos = 1 ; prot_pos <= nres ; ++prot_pos ) {
		ResidueType const & prot_type( pose.residue_type( prot_pos ) );
		if(  !prot_type.is_protein() ) continue;

		// map will automatically sort the "contacts" with the lowest total_score at the front of map
		std::map< Real, Size > contacts;

		// Loop over positions, skipping non-ligand
		for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
			ResidueType const & lig_type( pose.residue_type( lig_pos ) );
			if(  !lig_type.is_ligand() ) continue;

			Real pack_score = get_packing_score( pose, lig_pos, prot_pos );
			Real hb_score = get_hbond_score( pose, lig_pos, prot_pos );
			Real water_score = get_water_hbond_score( pose, lig_pos, prot_pos );

			Real total_score = pack_score + hb_score + water_score;
			if( pack_score < -0.5 || hb_score < -0.3 || water_score < -0.3 ) {
				contacts[total_score] = lig_pos;
				std::cout << "Energies between " << prot_pos << " and " << lig_pos << "are total: " << total_score << " pack: " << pack_score << " hbond: " << hb_score << " water: " << water_score << std::endl;
			}

		} // End loop over ligand residues

		Size contactssize( contacts.size() );
		if( contactssize != 0 ) {
			std::vector< Size > final_contacts;
			Size contactssize( contacts.size() );
			for( std::map< Real, Size >::const_iterator it( contacts.begin() ),
					end( contacts.end() ); it != end; ++it ) {
				if( contactssize == 1 ) {
					final_contacts.push_back( it->second );
					break;
				}
				Real first( (it)->first );
				Real firstb( (++it)->first );
				Real divided( first / firstb );
				if( divided > 1.5 ) {
					final_contacts.push_back( it->second );
					break;
				} else {
					final_contacts.push_back( it->second );
					final_contacts.push_back( (--it)->second );
					break;
				}
			} // End loop over the inital map of potential contacts, getting rid of weaker contacts

			utility::vector1< Size > final_final_contacts;

			int prot_pos_pdb( pose.pdb_info()->number( prot_pos ) );
			char prot_pos_chain( pose.pdb_info()->chain( prot_pos ) );
			std::string delimiter( "_" );

			for( Size ic=0; ic < final_contacts.size(); ++ic ) {
				if( pose.residue( prot_pos ).name3() == "GLY" ) continue;
				//protocols::motifs::Motif motif( pose.residue( prot_pos ), pose.residue( final_contacts[ic] ) );
				//conformation::ResidueOP protres = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( motif.restype_name1() ) );
				//conformation::ResidueOP dnares = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
				//conformation::ResidueOP dnares2 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
				//motif.place_residue( *protres, *dnares );
			//	bool broke( false );

				Size dna_pos = final_contacts[ ic ];
				int dna_pos_pdb( pose.pdb_info()->number( dna_pos ) );
				char dna_pos_chain( pose.pdb_info()->chain( dna_pos ) );
				std::string motif_name = pose.residue_type( prot_pos ).name3() + string_of( prot_pos_pdb ) + string_of( prot_pos_chain) + delimiter +  pose.residue_type( dna_pos ).name1() + string_of( dna_pos_pdb ) + string_of( dna_pos_chain ) + delimiter + pdb_name;
				//motif.store_remark( motif_name );

		//		for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motifs.begin(), end_itr = motifs.end();
		//				motifcop_itr != end_itr; ++motifcop_itr ) {
		//			protocols::motifs::MotifCOP motifcop( *motifcop_itr );
		//			if( motifcop->restype_name1() != motif.restype_name1() ) continue;
		//			if( motifcop->restype_name2() != motif.restype_name2() ) continue;
		//			motifcop->place_residue( *protres, *dnares2 );
		//			Real rmsdtest = scoring::automorphic_rmsd( *dnares, *dnares2, false );
					// will get appropriate name and add to motifs as well as outputting to file
					//std::cout << "RMSD = " << rmsdtest << " for motif at positions " << prot_pos << " and " << final_contacts[ic] << std::endl;
		//			if( rmsdtest < 0.2 ) {
		//				std::cout << "Skipping motif " << motif.remark() << " because it matches motif " << motifcop->remark() << " already in motif library, with an RMSD = " << rmsdtest << std::endl;
		//				broke = true;
			//			break;
			//		}
				}
				//if( !broke ) {
				//	std::cout << "Adding motif " << motif.remark() << std::endl;
				//	final_final_contacts.push_back( final_contacts[ic] );
			//		motifs.add_to_library( motif );
			//	}
			//}  // End loop that adds motifs to the MotifLibrary and skips the motif if a very similar motif has already been found
			//if( final_final_contacts.size() >= 1 ) {
			//	output_single_motif( pose, pdb_name, prot_pos, final_final_contacts );
			//}
		} // if there were actually contacts found for the particular protein residue
	} // End loop over protein residues

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
process_file_list()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace io::pdb;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace optimization;

  ScoreFunction scorefxn;

  scorefxn.set_weight( fa_atr, 1.00 );
  scorefxn.set_weight( fa_rep, 1.00 );
  scorefxn.set_weight( fa_sol, 1.00 );
  scorefxn.set_weight( fa_dun, 1.00 );
  scorefxn.set_weight( fa_pair, 1.00 );
  scorefxn.set_weight( p_aa_pp, 1.00 );
  scorefxn.set_weight( hbond_bb_sc, 1.0 );
  scorefxn.set_weight( hbond_sc, 1.0 );

	utility::vector1< std::string > pdb_files( start_files() );
	protocols::motifs::MotifLibrary motifs;
	for ( utility::vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
	      pdb_file != pdb_files.end(); ++pdb_file ) {
		std::string pdb_name( *pdb_file );

		std::cout << "Working on file: " << pdb_name << std::endl;

		std::string pdb_name4( pdb_name, pdb_name.size() - 8, 4 );

		Pose pose;
		core::import_pose::pose_from_file( pose, pdb_name , core::import_pose::PDB_file);

		process_for_motifs( pose, pdb_name4, motifs );
	}

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{

	try {

	devel::init( argc, argv );

	process_file_list();


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
