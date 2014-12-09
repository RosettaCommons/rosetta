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
//#include <core/options/option.hh>

#include <sstream>
#include <string>
#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>

////////////////////////////
////////////////////////////
////////////////////////////
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoreType.hh>

///////////////////////////////
//new, from upstream filter
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>
/////////////////////////////////////////

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
//#include <core/pack/types.hh>

#include <core/graph/Graph.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Map.Pose.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/scoring/mm/MMTorsionLibrary.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>

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
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/dna/util.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/Motif.hh>

#include <core/import_pose/import_pose.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

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
	protocols::motifs::MotifLibrary & motifs,
	std::string & pdb_name,
	Size prot_pos,
	Size lig_pos,
	utility::vector1< Size > & lig_atoms
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
		conformation::Residue protres( src_pose.residue( prot_pos ) );
		conformation::Residue ligres( src_pose.residue( lig_pos ) );
		//Here is where we deal with the 3 ligand motif atoms
		 for(core::Size atom_i = 1; atom_i <= ligres.natoms(); ++atom_i) {
			if(atom_i != lig_atoms[1]	&& atom_i != lig_atoms[2] && atom_i != lig_atoms[3]){
					numeric::xyzVector< core::Real > newpos(ligres.xyz(atom_i));
					newpos[1] = newpos[1]+100;
					ligres.set_xyz(atom_i, newpos);
				}
		 }
		Pose pose;
		pose.append_residue_by_jump( protres, 1 );
		pose.append_residue_by_jump( ligres, 1 );
		pose.conformation().insert_chain_ending( 1 );
		protocols::motifs::Motif motif(protres, ligres, lig_atoms);
		motifs.add_to_library( motif ); //Add motif to library
		int lig_pos_pdb( src_pose.pdb_info()->number( lig_pos ) );
		char lig_pos_chain( src_pose.pdb_info()->chain( lig_pos ) );
		std::string lig_talk = "Ligatoms" + delimiter + string_of( lig_atoms[1] )  + delimiter + string_of( lig_atoms[2] )  + delimiter + string_of( lig_atoms[3]) ;
		std::string motif_file_name_2 = motif_file_name + src_pose.residue_type( lig_pos ).name1() + string_of( lig_pos_pdb ) +  delimiter + lig_talk + delimiter +  string_of( lig_pos_chain ) + delimiter;

		std::string motif_file_name_3 = motif_file_name_2 + pdb_name + extension;

		std::cout << "Writing " << motif_file_name_3 << std::endl;
		io::pdb::dump_pdb( pose, motif_file_name_3 );
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

////////////////////////////////////////
//From protocols/match/output/UpstreamCollisionFilter.cc

    using namespace core::scoring;
    using namespace core::scoring::etable;
    using namespace core::scoring::methods;
    EnergyMethodOptions eopts;

    core::scoring::methods::ShortRangeTwoBodyEnergyOP etable_energy_ = new TableLookupEtableEnergy(
      *(ScoringManager::get_instance()->etable( eopts.etable_type() )), eopts );
   using namespace core::scoring;
    EnergyMap emap;
       // emap[ fa_atr ] = 0; emap[ fa_rep ] = 0; emap[ fa_sol ] = 0;
        etable_energy_->residue_pair_energy(
					pose.residue( pos1 ), pose.residue( pos2 ),
          pose, pack_score,
          emap );

//        Real energy = emap[ fa_atr ] + emap[ fa_rep ] + emap[ fa_sol ];
        Real energy = emap[ fa_atr ];

	return energy;
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
	std::cout << "in process_for_motifs function" << std::endl;

//////////////////////
//////////////////////
//////////////////////

		for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
			ResidueType const & lig_type( pose.residue_type( lig_pos ) );

// .is_ligand() comes from  src/core/chemical/ResidueType.cc (I think)
			if(  !lig_type.is_ligand() ) continue;
	    std::cout << "in ligand splitter block, found my ligand, lig_pos is " << lig_pos << std::endl;

 // This is to make a ligres object once we find our ligand
      conformation::Residue ligres( pose.residue( lig_pos ) );

			utility::vector1< utility::vector1< Size > > motif_indices_list;
			utility::vector1< utility::vector1< Size > > all_motif_indices;


        for(core::Size atom_i = 1; atom_i <= ligres.natoms(); ++atom_i) {
if( ligres.atom_type(atom_i).is_hydrogen()) {
						break;
																					}

	    std::cout << "in atom iterate block, atom num is" << atom_i <<   std::endl;
      std::string const atom_name = ligres.atom_name(atom_i);
	    std::cout << "atom name is " << atom_name <<  std::endl;

	// This is a for loop to iterate over each atom's connected atoms:
	core::conformation::Residue::AtomIndices atom_i_connects(  ligres.bonded_neighbor( atom_i ) );
 	for(core::Size atom_j = 1; atom_j <= atom_i_connects.size(); ++ atom_j ) {
if( ligres.atom_type(atom_i_connects[atom_j]).is_hydrogen()) {
						break;
																					}
		std::cout << "ATOM j: " << atom_i_connects[atom_j] << " Name: " << ligres.atom_name(atom_i_connects[atom_j]) << std::endl;

			// This is the next for loop to find connects for the second atom, giving us the final atom number (atom k)
			core::conformation::Residue::AtomIndices atom_j_connects(  ligres.bonded_neighbor( atom_i_connects[atom_j] ) );
	 		for(core::Size atom_k = 1; atom_k <= atom_j_connects.size(); ++ atom_k ) {
if( ligres.atom_type(atom_j_connects[atom_k]).is_hydrogen()) {
						break;
																															}
		std::cout << "ATOM k: " << atom_j_connects[atom_k] << " Name: " << ligres.atom_name(atom_j_connects[atom_k]) << std::endl;

				chemical::AtomType atom_i_type(ligres.atom_type(atom_i));
				std::string atom_i_name = atom_i_type.atom_type_name();

  		   std::cout << "Connected triplet is: " << atom_i << ", type is " << atom_i_name  << ", ";
    		 std::cout << atom_i_connects[atom_j] << ", type is " << ligres.atom_type_index(atom_i_connects[atom_j]) << ", " ;
     	   std::cout << atom_j_connects[atom_k] << ", type is " << ligres.atom_type_index(atom_j_connects[atom_k]) << " " << std::endl;
				 if ( atom_i != atom_j_connects[atom_k] ) {

						//make the 3 atom vector
							utility::vector1< Size > cur_motif_indices;
							cur_motif_indices.push_back( atom_i );
							cur_motif_indices.push_back( atom_i_connects[atom_j] );
							cur_motif_indices.push_back( atom_j_connects[atom_k] );

						//check if current index is a duplicate

						bool current_is_duplicate ( false );
						if ( !motif_indices_list.empty() ) {
							for(core::Size cur_motif_check = 1; cur_motif_check <= motif_indices_list.size(); ++ cur_motif_check) {
								// run a test to see if vectors are identical
								utility::vector1< Size > cur_mainlist ( motif_indices_list[cur_motif_check] );
								utility::vector1< Size > sort_from_mainlist (cur_mainlist) ;
								    std::sort(  sort_from_mainlist.begin(),  sort_from_mainlist.end() );
								utility::vector1< Size > sort_from_curindex (cur_motif_indices) ;
								    std::sort( sort_from_curindex.begin(), sort_from_curindex.end() );
								if ( sort_from_mainlist[1] == sort_from_curindex[1] && sort_from_mainlist[2] == sort_from_curindex[2] && sort_from_mainlist[3] == sort_from_curindex[3] ) {
									current_is_duplicate = true;
																																												 }
																																											                              }
																								 }
						//Making master list to see what's getting pruned with my check
									 all_motif_indices.push_back( cur_motif_indices );
						//if current index isn't a duplicate, add it to the vector
								if ( !current_is_duplicate ) {
									 motif_indices_list.push_back( cur_motif_indices );
																				 }
                                 }

}

}
   std::cout << std::endl;

				}


 std::cout << "Total 3 atoms in pruned indices list is: " << motif_indices_list.size()  << std::endl;
 std::cout << "Total 3 atoms in unpruned indices list is: " << all_motif_indices.size()  << std::endl;

////////////////////////////
////////////////////////////
////////////////////////////


	for( int prot_pos = 1 ; prot_pos <= nres ; ++prot_pos ) {
		ResidueType const & prot_type( pose.residue_type( prot_pos ) );
		if(  !prot_type.is_protein() ) continue;
		if( pose.residue( prot_pos ).name3() == "GLY" ) continue;

		// map will automatically sort the "contacts" with the lowest total_score at the front of map
		std::map< Real, Size > contacts;
		std::map< Real, Size  > distance_sorter;

		// Loop over positions, skipping non-ligand
		for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
			ResidueType const & lig_type( pose.residue_type( lig_pos ) );

			if(  !lig_type.is_ligand() ) continue;

			Real pack_score = get_packing_score( pose, lig_pos, prot_pos );

			if( pack_score < -0.5 ) { //Enter into motif checking loop--For each interacting residue:
				contacts[pack_score] = lig_pos;
				std::cout << "Residue " << prot_pos << " passed energy cut with " << pack_score << std::endl;

				for(core::Size motif_position = 1; motif_position <= motif_indices_list.size(); ++ motif_position) {	//for each 3 atom triplet from ligand
				 bool resi_trip_match( false );
				 utility::vector1< Size > cur_trip ( motif_indices_list[motif_position] );
				 Real closest_distance(5.0);
				 for(core::Size cur_trip_pos = 1; cur_trip_pos <= cur_trip.size(); ++ cur_trip_pos) {  	//for each atom in ligand triplet
					for(core::Size residue_atom_number = 1; residue_atom_number <= pose.residue( prot_pos ).nheavyatoms(); ++ residue_atom_number) {   //for each atom in current residue
					  Real atom_atom_distance( pose.residue( lig_pos ).xyz( cur_trip[cur_trip_pos]  ).distance( pose.residue( prot_pos ).xyz( residue_atom_number ) ) );
						if( atom_atom_distance < 4.0 ) {
							resi_trip_match = true;
							if( atom_atom_distance < closest_distance) { closest_distance=atom_atom_distance; }
						}
		}
		}
			//here is the end of what we do with each triplet, this is where we add triplet-residue to vector of pairs
			if( resi_trip_match ) {
							distance_sorter[closest_distance] = motif_position;
      }
		}
//Here is where we find the closest few triplets for the current residue (we are in a protein residue loop)

      utility::vector1< Size > top_triplets;
    Size sorter_size( distance_sorter.size() );

/////////////////////
/////////////////////
    if( sorter_size != 0 ) {
					int min(15); //This is the max number of motifs to output
// std::cout << "Minimum is " << min << std::endl;
	Size stop_after_min_count(1);
////////////////////////////////////////////////////////////////////////////////////////////
      for( std::map< Real, Size >::const_iterator it( distance_sorter.begin() ),  end( distance_sorter.end() ); it != end; ++it ) {
					if( stop_after_min_count == min) { break; }
					Size my_cur_trip( it-> second );
					top_triplets.push_back( my_cur_trip );
					++stop_after_min_count;
				}
///////////////////////////////////////////////////////////////////////////////////////
        }
///////////////////
///////////////////

        std::cout << "Top triplets contains " << top_triplets.size() << " items." << std::endl << "Top triplets are: " ;
 for(core::Size top_trip_pos = 1; top_trip_pos <= top_triplets.size(); ++ top_trip_pos) {
 					Size this_triplet_number(top_triplets[top_trip_pos]);
					utility::vector1< Size > this_triplet( motif_indices_list[this_triplet_number] );
					std::cout <<  this_triplet_number << ": " << this_triplet[1] << "-" << this_triplet[2] << "-" <<  this_triplet[3]  ;
					output_single_motif( pose, motifs, pdb_name, prot_pos, lig_pos, this_triplet ); // Output a pdb for each protres -- ligand motif pair.  This function also makes motif for current triplet/resi pair and adds it to the motif library.
			}
std::cout <<  std::endl;
std::cout <<  std::endl;
		} // End "if good energy" statement

/*
		Size contactssize( contacts.size() );
		if( contactssize != 0 ) {
			std::vector< Size > final_contacts; //NEVER DO THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
		} // if there were actually contacts found for the particular protein residue */ //Here is the end of commenting out the useless section
	} // End loop over protein residues
}
//Here we're going to check to see what's in motif_indices_list
for(core::Size motif_position = 1; motif_position <= motif_indices_list.size(); ++ motif_position) {
utility::vector1< Size > cur_trip ( motif_indices_list[motif_position] );
std::cout << "Motif index contains: " << cur_trip[1] << "-" << cur_trip[2] << "-" << cur_trip[3] << std::endl;
}
}
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

	//Here is where I create the MotifLibrary object which I populate with the motifs
	protocols::motifs::MotifLibrary motifs;
	for ( utility::vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
	      pdb_file != pdb_files.end(); ++pdb_file ) {
		std::string pdb_name( *pdb_file );

		std::cout << "Working on file: " << pdb_name << std::endl;

		//std::string pdb_name4( pdb_name, pdb_name.size() - 8, 4 );

		Pose pose;
   try{
		core::import_pose::pose_from_pdb( pose, pdb_name );
				    } catch( utility::excn::EXCN_BadInput excn ) {
											std::cout << "Got stuck in missing heavyatom loop, continuing" <<  std::endl;
											continue; }

		process_for_motifs( pose, pdb_name, motifs );

	} //End of process_for_motifs loop

	  protocols::motifs::MotifCOPs motifcops = motifs.library();

					  std::string filename( "AllMattMotifsFile.motifs" );
						  utility::io::ozstream motif_output_file( filename );
							  for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motifcops.begin(), end_itr = motifcops.end();
								      motifcop_itr != end_itr; ++motifcop_itr ) {
											      protocols::motifs::MotifCOP motifcop( *motifcop_itr );
														      motif_output_file << *motifcop;
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
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
