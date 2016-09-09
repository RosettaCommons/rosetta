// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file hydrophobic_patch_stats_calculator.cc
/// @brief A simple protocol which outputs information about surface residues on the passed in list of PDBs.
/// @author Ron Jacak

// Unit headers - These are headers which declare the things this file defines. In other words, A.cc includes A.hh here.

// Package Headers - These are closely related headers, within a subdirectory of core at least.
//                   For example, a packer-related header might have the PackerTask here.

// Project Headers - This is where loose headers from all over the core project go (scorefunction, packer, etc).
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/sasa.hh>

// Utility headers - Put things from the utility library here, as well as tracers and the option system.
#include <devel/init.hh>
#include <basic/tracer.hh>
#include <basic/database/open.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/viewer/viewers.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

#include <core/types.hh>
// Vector is a typedef defined in core/types.hh as follows: numeric::xyzVector< Length >  Vector;
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>

#include <utility/exit.hh>
#include <basic/prof.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <sstream>
#include <string>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>

//#include <math>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ronj.hydrophobic_patch_stats_calculator" );

using utility::file::FileName;
using utility::file::file_basename;
using utility::file::trytry_ofstream_open;
using utility::file::file_exists;

using namespace core;
using namespace ObjexxFCL::format;

// make the "surface-exposed" cut-off a variable
core::Size SURFACE_EXPOSED_CUTOFF = 20;

// any residue with more than 24 nbs will automatically get 0.00 for the hASA
core::Size BURIED_RESIDUE_NO_HSASA_CUTOFF = 24;

#define FILE_DEBUG 1

/// @brief Calculate the weighted neighbor count given an upper and lower bound

/// @details
// Copied straight from NVscore.cc written by Sam DeLuca
//
core::Real
neighborWeight( core::Vector::Value& dist, core::Real lBound, core::Real uBound ) {

	if ( dist <= lBound ) {
		return 1;
	} else if ( (lBound < dist) && (uBound > dist) ) {
		//if between upper and lower bound, score follows a smooth function
		return ( cos( ( (dist-lBound) / (uBound-lBound) ) * numeric::constants::r::pi ) + 1 )/2.0;
	}

	return 0;  // else dist >= uBound
}


/// @details
/// Read in the values contained in the NACCESS data file into the map passed in by reference.
///
void
read_NACCESS_asa_file( std::string filename, std::map< std::string, Real > & map ) {

	std::ifstream naccess_data_ifstream( filename.c_str() );
	while ( naccess_data_ifstream.good() ) {
		std::string line, key;
		Real asa = 0.0;

		while ( !naccess_data_ifstream.eof() ) {
			getline( naccess_data_ifstream, line );
			std::istringstream linestream( line );
			linestream >> key;
			linestream >> asa;
			map[ key ] = asa;
		}
	}

}


/// @brief a quick protocol written to calculate some statistics on surface residues in native and designed proteins
void
find_hppatches_nb_graph( std::vector< FileName > & pdb_file_names ) {

	// keep track of the number of times a patch of size [array_index] is seen
	std::vector<int > patch_size_counts_nb_graph;
	patch_size_counts_nb_graph = std::vector<int>(11,0);

	int totalNumPatches = 0;
	std::vector<std::pair<std::string, int> > countSEHydrophobic( 0 );

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		//TR << "Read in pdb file '" << pose.name() << "'" << std::endl;

		// create a score function using the standard packer weights
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

		scorefxn->score( pose );
		//TR << pose.energies() << std::endl;

		int patchSize = 0;
		int countSEHydrophobicPerStructure = 0;

		// To determine how often a surface hydrophobic neighbors other hydrophobics, we first have to loop over all
		// residues and find exposed hydrophobics
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		assert ( pose.size() > 4 );

		countSEHydrophobicPerStructure = 0;

		for ( core::Size res1_position = 4; res1_position <= pose.size() - 3; ++res1_position ) {

			// reset the counters, before EVERY iteration
			patchSize = 0;

			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

			// our definition of surface residue is that the residue has fewer than 16 neighbors
			if ( (Size)tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
				continue;
			}

			// passed the surface-exposed check...
			// check if this residue is a hydrophobic one. we want to consider all surface residues, not just hydrophobic ones
			// if it is, increment the hp count since we want to count all the hydrophobics

			// the is_polar() member function of type Residue return true for the following residue types:
			// ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR
			// the negative of is_polar() would return true for the following residues:
			// ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, VAL, TRP, TYR
			// reordered as below is
			// VAL, ILE, LEU, MET, PHE, GLY, ALA, PRO, TRP, TYR, and CYS
			// the residues I want to include as hydrophobic include:
			// VAL, ILE, LEU, MET, PHE, GLY, ALA, PRO, TRP, TYR
			// so, only CYS is included as a hydrophobic when it shouldn't be

			if ( pose.residue( res1_position ).is_polar() ) {
				continue;
			} else if ( !(pose.residue( res1_position ).is_polar()) && pose.residue( res1_position ).name3() != "CYS" ) {
				countSEHydrophobicPerStructure++;
				patchSize++;
			} else {  // it must be a "hydrophobic" CYS, so just go on
				continue;
			}

			//TR << "Neighbors of residue " << pose.residue( res1_position ).name3() << " " << res1_position << " include " << std::endl;
			// for every Edge in the neighbor graph, figure out if that residue is surface exposed *and* hydrophobic
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_end(); eli != eli_end; ++eli ) {

				// save the value to simplify code ahead
				int res2_position = (*eli)->get_other_ind( res1_position );

				// get the other node for this edge, so pass in the res1 node to this method
				//TR << pose.residue( res2_position ).name3() << " " << res2_position << std::endl;

				if ( !(pose.residue( res2_position ).is_polar()) && (pose.residue( res2_position ).name3() != "CYS") ) {
					// ok, so it's hydrophobic and not CYS, but is it surface-exposed, too?
					if ( (Size)tenA_neighbor_graph.get_node( res2_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
						continue;
					}
					// passed all checks
					patchSize++;
				}
			}

			// now that we know how many surface-exposed, hphobic neighbors res1 has, save it in an array
			if ( patchSize != 0 ) {
				totalNumPatches++;
				if ( patchSize > 10 ) {
					patchSize = 10;  // don't keep stats above 10 at this point
				}
				patch_size_counts_nb_graph[patchSize] = patch_size_counts_nb_graph[patchSize] + 1;
			}

		} // end res1 loop

		// store the value outside the threshold loop so we don't get 5 repeated lines
		countSEHydrophobic.push_back( std::pair<std::string,int>( pose.pdb_info()->name(), countSEHydrophobicPerStructure ) );

		// output pdb file
		//core::io::pdb::dump_pdb( pose, file_basename( pose.name()) + "_hp.pdb" );

	} // end for loop over multiple input pdb files

	TR << "Total patches found using neighbor graph method: " << std::endl;
	TR << patch_size_counts_nb_graph.size() << std::endl;

} // end find_hppatches_nb_graph


/// @brief a quick protocol written to calculate some statistics on surface residues in native and designed proteins
void
find_hppatches_distance( std::vector< FileName > & pdb_file_names ) {

	// keep track of the number of times a patch of size [array_index] is seen
	// so patch_size_counts[6][4]=23 means a patch of size 4 was seen 23 times in the protein using a
	// threshold radius of 6A.
	//std::map<int, std::vector<int > > patch_size_counts( (6, std::vector<int>(11, 0)), (7, std::vector<int>(11, 0)),
	//				(8, std::vector<int>(11, 0)), (9, std::vector<int>(11, 0)), (10, std::vector<int>(11, 0)) );
	std::map<int, std::vector<int > > patch_size_counts;
	patch_size_counts[6] = std::vector<int>(11,0);
	patch_size_counts[7] = std::vector<int>(11,0);
	patch_size_counts[8] = std::vector<int>(11,0);
	patch_size_counts[9] = std::vector<int>(11,0);
	patch_size_counts[10] = std::vector<int>(11,0);

	std::map<int,int> totalNumPatches; //(6,0), (7,0), (8,0), (9,0), (10,0)
	totalNumPatches[6] = 0;
	totalNumPatches[7] = 0;
	totalNumPatches[8] = 0;
	totalNumPatches[9] = 0;
	totalNumPatches[10] = 0;

	std::vector<std::pair<std::string, int> > countSEHydrophobic( 0 );

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		//TR << "Read in pdb file '" << pose.name() << "'" << std::endl;

		// create a score function using the standard packer weights
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

		scorefxn->score( pose );
		//TR << pose.energies() << std::endl;


		int patchSize = 0;
		int countSEHydrophobicPerStructure = 0;

		// To determine how often a surface hydrophobic neighbors other hydrophobics, we first have to loop over all
		// residues and find exposed hydrophobics
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		assert ( pose.size() > 4 );

		int thresholdRadius;
		for ( thresholdRadius = 6; thresholdRadius <= 10; thresholdRadius++ ) {

			countSEHydrophobicPerStructure = 0;

			for ( core::Size res1_position = 4; res1_position <= pose.size() - 3; ++res1_position ) {

				// reset the counters, before EVERY iteration
				patchSize = 0;

				core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

				// our definition of surface residue is that the residue has fewer than 16 neighbors
				if ( (Size)tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
					continue;
				}

				// passed the surface-exposed check...
				// check if this residue is a hydrophobic one. we want to consider all surface residues, not just hydrophobic ones
				// if it is, increment the hp count since we want to count all the hydrophobics

				// the is_polar() member function of type Residue return true for the following residue types:
				// ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR
				// the negative of is_polar() would return true for the following residues:
				// ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, VAL, TRP, TYR
				// reordered as below is
				// VAL, ILE, LEU, MET, PHE, GLY, ALA, PRO, TRP, TYR, and CYS
				// the residues I want to include as hydrophobic include:
				// VAL, ILE, LEU, MET, PHE, GLY, ALA, PRO, TRP, TYR
				// so, only CYS is included as a hydrophobic when it shouldn't be

				if ( pose.residue( res1_position ).is_polar() ) {
					continue;
				} else if ( !(pose.residue( res1_position ).is_polar()) && pose.residue( res1_position ).name3() != "CYS" ) {
					countSEHydrophobicPerStructure++;
					patchSize++;
				} else {  // it must be a "hydrophobic" CYS, so just go on
					continue;
				}


				// now find all other surface-exposed hydrophobic residues.
				for ( core::Size res2_position = 4; res2_position <= pose.size() - 3; ++res2_position ) {
					// see note above for starting/terminating conditions

					if ( res1_position == res2_position )
						continue;

					if ( (Size)tenA_neighbor_graph.get_node( res2_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF )
						continue;

					// if it's a polar residues, move on
					// should we include GLY as a polar residue, or not?  we'll see how that affects the stats later
					if ( pose.residue( res2_position ).is_polar() )
						continue;

					// we need to determine whether the centroid of the side chain (is that the same as the action coordinate?)
					// of this residue has a distance less than the centroid of the res1 sidechain.  This method is a brute
					// force way of determining "neighbor".  Another way would be to use the neighbor graph.  Try both and
					// make sure the results are the same.
					conformation::Residue const & rsd1 = pose.residue( res1_position );
					conformation::Residue const & rsd2 = pose.residue( res2_position );

					//float distanceBetweenAtoms =  rsd1.xyz(rsd1.nbr_atom()).distance( rsd2.xyz(rsd2.nbr_atom()) );
					//TR << "Distance between atoms " << rsd1.name3() << " " << res1_position
					//	<< " " << rsd1.atom_name( rsd1.nbr_atom() ) << " - " << rsd2.atom_name( rsd2.nbr_atom() )
					//	<< rsd2.name3() << " " << res2_position << ": " << distanceBetweenAtoms << std::endl;

					if ( rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz(rsd2.nbr_atom()) ) < thresholdRadius ) {

						// increment the res1 neighbor count variable, to be printed later
						// but we want to bin the neighbor counts depending on how surface-exposed the residue is
						patchSize++;
					}

				} // end res2 loop

				// now that we know how many surface-exposed, hphobic neighbors res1 has, save it in an array
				// in the pdb_stats_hppatch namespace, a 2D vector is defined and externed
				if ( patchSize != 0 ) {
					totalNumPatches[thresholdRadius]++;
					if ( patchSize > 10 ) {
						patchSize = 10;  // don't keep stats above 10 at this point
					}
					patch_size_counts[thresholdRadius][patchSize] = patch_size_counts[thresholdRadius][patchSize] + 1;
				}

			} // end res1 loop

			// store the value outside the threshold loop so we don't get 5 repeated lines
			countSEHydrophobic.push_back( std::pair<std::string,int>( pose.pdb_info()->name(), countSEHydrophobicPerStructure ) );

			// output pdb file
			//core::io::pdb::dump_pdb( pose, file_basename( pose.name()) + "_hp.pdb" );

		} // end thresholdRadius

	}  // end for loop over multiple input pdb files

	TR << "Total patches found using distance method: " << std::endl;
	TR << patch_size_counts.size() << std::endl;

} // end find_hppatches_distance


/// @brief a quick protocol written to calculate some statistics on surface residues in native and designed proteins
void
calculate_percent_hydrophobic_stats( std::vector< FileName > & pdb_file_names )
{

	// create a score function using the standard packer weights
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of percent hydrophobic for a given number of surace neighbors
	std::map<int, std::vector<float> > percent_hp_in_se_nbs_by_nb;
	std::map<int, std::vector<float> > percent_hp_in_all_nbs_by_nb;

	TR << "Processing... " << std::endl;
	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		TR << pose.pdb_info()->name() << std::endl;

		scorefxn->score( pose );

		// To determine how often a surface hydrophobic neighbors other hydrophobics, we first have to loop over all
		// residues and find exposed hydrophobics
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		assert ( pose.size() > 4 );
		for ( core::Size res1_position = 4; res1_position <= pose.size() - 3; ++res1_position ) {

			// these variables are per res1 values, reset every new exposed res1
			int countNeighbors = 0;
			int countSurfaceNeighbors = 0;
			int countHydrophobicSurfaceNeighbors = 0;

			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

			// our definition of surface residue is that the residue has fewer than 16 neighbors
			if ( (Size)tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// passed the surface-exposed check...
			// save the number of neighbors
			countNeighbors = tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() - 1;

			// now determine the number of *surface* neighbors
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_end(); eli != eli_end; ++eli ) {

				// save the value to simplify code ahead
				int res2_position = (*eli)->get_other_ind( res1_position );

				if ( (Size)tenA_neighbor_graph.get_node( res2_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
					// it's not surface-exposed, move on to the next one
					continue;
				} else {
					// it is surface-exposed
					countSurfaceNeighbors++;

					// now, also check if it's hphobic
					if ( !(pose.residue( res2_position ).is_polar()) && (pose.residue( res2_position ).name3() != "CYS") ) {
						// it is hydrophobic and not CYS, cool
						countHydrophobicSurfaceNeighbors++;
					}
				}

			} // end going through all neighbor graph edges

			//TR << "Residue " << pose.residue( res1_position ).name3() << " " << res1_position
			//	<< " is surface-exposed and has " << countNeighbors << " neighbors, " << countSurfaceNeighbors << " surface neighbors, and " <<
			//	countHydrophobicSurfaceNeighbors << " hphobic surface neighbors. " << std::endl;

			if ( countSurfaceNeighbors == 0 ) {
				// res1 has nb's, but no surface neighbors. must be kinda buried. don't tally it's count then.
				continue;
			}

			// now we've set the count of neighbors, surface nb's, and hphobic surface nb's.
			float percent_hphobic = ((float)countHydrophobicSurfaceNeighbors / countSurfaceNeighbors) * 100;
			percent_hp_in_se_nbs_by_nb[ countSurfaceNeighbors ].push_back( percent_hphobic );

			percent_hphobic = ((float)countHydrophobicSurfaceNeighbors / countNeighbors) * 100;
			percent_hp_in_all_nbs_by_nb[ countSurfaceNeighbors ].push_back( percent_hphobic );

		} // end res1 loop

	} // end for loop over multiple input pdb files

	TR << "Stats calculated: " << std::endl;

	// use const_iterator to walk through elements of pairs
	for ( std::map<int, std::vector<float> >::const_iterator iter = percent_hp_in_se_nbs_by_nb.begin(); iter != percent_hp_in_se_nbs_by_nb.end(); ++iter ) {
		//std::cout << iter->first << '\t' << iter->second << std::endl;
	}

	// average
	std::map<int, float> average_percent_hp_in_se_nbs_by_nb;
	for ( std::map<int, std::vector<float> >::const_iterator iter = percent_hp_in_se_nbs_by_nb.begin(); iter != percent_hp_in_se_nbs_by_nb.end(); ++iter ) {
		float sum = 0;
		for ( Size i=0; i < iter->second.size(); ++i ) {
			sum = sum + iter->second[i];
		}
		float average = sum / iter->second.size();
		average_percent_hp_in_se_nbs_by_nb[ iter->first ] = average;
	}
	TR << "average_percent_hp_in_se_nbs_by_nb" << std::endl;
	for ( std::map<int, float>::const_iterator iter = average_percent_hp_in_se_nbs_by_nb.begin(); iter != average_percent_hp_in_se_nbs_by_nb.end(); ++iter ) {
		//std::cout << iter->first << '\t' << iter->second << std::endl;
	}

	// standard deviation
	std::map<int, float> std_dev_hp_in_se_nbs_by_nb;
	for ( std::map<int, std::vector<float> >::const_iterator iter = percent_hp_in_se_nbs_by_nb.begin(); iter != percent_hp_in_se_nbs_by_nb.end(); ++iter ) {
		float sum_of_diffs = 0;
		for ( Size i=0; i < iter->second.size(); ++i ) {
			float diff = iter->second[i] - average_percent_hp_in_se_nbs_by_nb[ iter->first ];
			float diff_squared = diff * diff;
			sum_of_diffs = sum_of_diffs + diff_squared;
		}
		float std_dev = sqrt( sum_of_diffs / (iter->second.size()-1) );
		std_dev_hp_in_se_nbs_by_nb[ iter->first ] = std_dev;
	}
	TR << "std_dev_hp_in_se_nbs_by_nb" << std::endl;
	for ( std::map<int, float>::const_iterator iter = std_dev_hp_in_se_nbs_by_nb.begin(); iter != std_dev_hp_in_se_nbs_by_nb.end(); ++iter ) {
		std::cout << iter->second << std::endl;
	}

	// average
	std::map<int, float> average_percent_hp_in_all_nbs_by_nb;
	for ( std::map<int, std::vector<float> >::const_iterator iter = percent_hp_in_all_nbs_by_nb.begin(); iter != percent_hp_in_all_nbs_by_nb.end(); ++iter ) {
		float sum = 0;
		for ( Size i=0; i < iter->second.size(); ++i ) {
			sum = sum + iter->second[i];
		}
		float average = sum / iter->second.size();
		average_percent_hp_in_all_nbs_by_nb[ iter->first ] = average;
	}
	TR << "average_percent_hp_in_all_nbs_by_nb" << std::endl;
	for ( std::map<int, float>::const_iterator iter = average_percent_hp_in_all_nbs_by_nb.begin(); iter != average_percent_hp_in_all_nbs_by_nb.end(); ++iter ) {
		std::cout << iter->first << '\t' << iter->second << std::endl;
	}

	// standard deviation
	std::map<int, float> std_dev_hp_in_all_nbs_by_nb;
	for ( std::map<int, std::vector<float> >::const_iterator iter = percent_hp_in_all_nbs_by_nb.begin(); iter != percent_hp_in_all_nbs_by_nb.end(); ++iter ) {
		float sum_of_diffs = 0;
		for ( Size i=0; i < iter->second.size(); ++i ) {
			float diff = iter->second[i] - average_percent_hp_in_all_nbs_by_nb[ iter->first ];
			float diff_squared = diff * diff;
			sum_of_diffs = sum_of_diffs + diff_squared;
		}
		float std_dev = sqrt( sum_of_diffs / (iter->second.size()-1) );
		std_dev_hp_in_all_nbs_by_nb[ iter->first ] = std_dev;
	}
	TR << "std_dev_hp_in_all_nbs_by_nb" << std::endl;
	for ( std::map<int, float>::const_iterator iter = std_dev_hp_in_all_nbs_by_nb.begin(); iter != std_dev_hp_in_all_nbs_by_nb.end(); ++iter ) {
		std::cout << iter->second << std::endl;
	}


} // end calculate_percent_hydrophobic_stats


/// @brief a quick protocol written to calculate some statistics on surface residues in native and designed proteins
void
calculate_percent_hydrophobic_distribution( std::vector< FileName > & pdb_file_names )
{

	// create a score function using the standard packer weights
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of percent hydrophobic for a given number of surace neighbors
	std::vector<float> percent_hp_in_se_nbs;
	std::vector<int> count_hp_in_se_nbs;

	TR << "Processing... " << std::endl;
	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		TR << pose.pdb_info()->name() << std::endl;

		scorefxn->score( pose );

		// To determine how often a surface hydrophobic neighbors other hydrophobics, we first have to loop over all
		// residues and find exposed hydrophobics
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		assert ( pose.size() > 4 );
		for ( core::Size res1_position = 4; res1_position <= pose.size() - 3; ++res1_position ) {

			// these variables are per res1 values, reset every new exposed res1
			int countNeighbors = 0;
			int countSurfaceNeighbors = 0;
			int countHydrophobicSurfaceNeighbors = 0;

			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

			// our definition of surface residue is that the residue has fewer than 16 neighbors
			if ( (Size)tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// passed the surface-exposed check...
			// save the number of neighbors
			countNeighbors = tenA_neighbor_graph.get_node( res1_position )->num_neighbors_counting_self() - 1;

			// now determine the number of *surface* neighbors
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node( res1_position )->const_edge_list_end(); eli != eli_end; ++eli ) {

				// save the value to simplify code ahead
				int res2_position = (*eli)->get_other_ind( res1_position );

				if ( (Size)tenA_neighbor_graph.get_node( res2_position )->num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
					// it's not surface-exposed, move on to the next one
					continue;
				} else {
					// it is surface-exposed
					countSurfaceNeighbors++;

					// now, also check if it's hphobic
					if ( !(pose.residue( res2_position ).is_polar()) && (pose.residue( res2_position ).name3() != "CYS") ) {
						// it is hydrophobic and not CYS, cool
						countHydrophobicSurfaceNeighbors++;
					}
				}

			} // end going through all neighbor graph edges

			//TR << "Residue " << pose.residue( res1_position ).name3() << " " << res1_position
			//	<< " is surface-exposed and has " << countNeighbors << " neighbors, " << countSurfaceNeighbors << " surface neighbors, and " <<
			//	countHydrophobicSurfaceNeighbors << " hphobic surface neighbors. " << std::endl;

			if ( countSurfaceNeighbors == 0 ) {
				// res1 has nb's, but no surface neighbors. must be kinda buried. don't tally it's count then.
				continue;
			}

			// now we've set the count of neighbors, surface nb's, and hphobic surface nb's.
			float percent_hphobic = ((float)countHydrophobicSurfaceNeighbors / countSurfaceNeighbors) * 100;
			percent_hp_in_se_nbs.push_back( percent_hphobic );

			count_hp_in_se_nbs.push_back( countHydrophobicSurfaceNeighbors );

		} // end res1 loop

	} // end for loop over multiple input pdb files

	TR << "Stats calculated: " << std::endl;

	std::ofstream percent_file;
	if ( trytry_ofstream_open( percent_file, "percent.txt", std::ios::out ) ) {
		TR << "outputting percent_hp_in_se_nbs" << std::endl;
		for ( Size i = 0; i < percent_hp_in_se_nbs.size(); ++i ) {
			percent_file << percent_hp_in_se_nbs[i] << std::endl;
		}
	} else {
		TR << "Could not open filename 'percent.txt'." << std::endl;
	}

	std::ofstream counts_file;
	if ( trytry_ofstream_open( counts_file, "counts.txt", std::ios::out ) ) {
		TR << "outputting count_hp_in_se_nbs" << std::endl;
		for ( Size i = 0; i < count_hp_in_se_nbs.size(); ++i ) {
			counts_file << count_hp_in_se_nbs[i] << std::endl;
		}
	} else {
		TR << "Could not open filename 'counts.txt'." << std::endl;
	}


} // end calculate_percent_hydrophobic_distribution


/// @brief a quick protocol written to calculate hydrophobic ASA on a set of proteins

/// @details
// NACCESS generates 3 files for every PDB: a log file, an atomic accessibility file (.asa) and a residue accessibility file (.rsa)
// Instead of using the complicated C++ string parsing routines, I'll use filesystem grep calls to make a residue-level file
// which has the residue string/key as the first column and the total hydrophobic ASA for that residue.

// The first two residues from RalA structure, from the .asa file
// ATOM      1  N   SER A  11      84.265  65.648  63.661  41.833  1.65
// ATOM      2  CA  SER A  11      83.271  66.744  63.480  18.110  1.87
// ATOM      3  C   SER A  11      83.851  67.907  62.673   2.018  1.76
// ATOM      4  O   SER A  11      83.138  68.855  62.329  28.553  1.40
// ATOM      5  CB  SER A  11      82.779  67.237  64.848  46.516  1.87
// ATOM      6  OG  SER A  11      83.864  67.568  65.700  25.010  1.40
//
// ATOM      7  N   LEU A  12      85.145  67.835  62.374   1.075  1.65
// ATOM      8  CA  LEU A  12      85.801  68.876  61.590   5.463  1.87
// ATOM      9  C   LEU A  12      85.558  68.630  60.108   0.810  1.76
// ATOM     10  O   LEU A  12      85.292  67.499  59.697  16.143  1.40
// ATOM     11  CB  LEU A  12      87.311  68.878  61.849   1.968  1.87
// ATOM     12  CG  LEU A  12      87.844  69.433  63.172   6.956  1.87
// ATOM     13  CD1 LEU A  12      87.256  68.660  64.346  38.857  1.87
// ATOM     14  CD2 LEU A  12      89.365  69.337  63.171   1.469  1.87

// REM  Relative accessibilites read from external file "standard.data"
// REM  File of summed (Sum) and % (per.) accessibilities for
// REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
// REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
// RES SER A  11   162.04 139.1  89.64 114.8  72.40 188.6  66.64 137.3  95.40 140.4
// RES LEU A  12    72.74  40.7  54.71  38.8  18.03  48.1  55.52  39.0  17.22  47.4

// If we take the sum of the C(arbon) and S(ulfur) atoms in the .asa file, we get the total nonpolar (hp) ASA for that residue
// as seen in the .rsa file.  But that sum includes the carbonyl carbon in the amino acid.  Should that really contribute
// to the hydrophobic ASA of that residue?  It is true that every amino acid regardless of its side chain will have that
// carbonyl carbon so it shouldn't really have a discriminatory effect when placing different residues. The carbonyl carbon
// should be buried/exposed to the same degree regardless of residue type. The only thing I need to be careful about is when
// I go about calculating the total hydrophobic ASA for a protein that I also include this carbonyl carbon in the calculation.
// Since I'm going to be optimizing the score function to have a total hASA like that seen in natives using the average hASA
// of a residue, I need to be consistent and use the same approach for both.

void
calculate_hydrophobic_accessible_surface_area( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has, broken up by res type
	utility::vector1< utility::vector1< Real > > hp_ASA_by_aa_type_nbs_1_to_10( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hp_ASA_by_aa_type_nbs_11_to_13( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hp_ASA_by_aa_type_nbs_14_to_16( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hp_ASA_by_aa_type_nbs_17_to_20( chemical::num_canonical_aas );

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}


		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		scorefxn->score( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hp_ASA;

		std::ifstream naccess_data_ifstream( naccess_asa_filename.str().c_str() );
		while ( naccess_data_ifstream.good() ) {

			std::string line;
			std::string key;
			Real asa = 0.0;

			while ( !naccess_data_ifstream.eof() ) {
				getline( naccess_data_ifstream, line );

				std::istringstream linestream( line );
				linestream >> key;
				linestream >> asa;

				res_to_calculated_hp_ASA[ key ] = asa;
			}

		}

		// First, we have to figure out which residues are surface-exposed
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) continue;

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			Size countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// passed the surface-exposed check...
			// something else to be careful about is that mini likes to renumber the residues that it reads from a PDB file
			// so if the PDB starts at something other than 1, we need to translate the mini number to the PDB number

			// convert the PDB resid to the pose resid in case the PDB numbering is screwy
			int pdb_resid = pose.pdb_info()->number(ii);
			char icode = pose.pdb_info()->icode(ii);

			//if ( pdb_resid != ii ) {
			//	TR << "Converted mini res num " << ii << " to PDB res num " << pdb_resid << std::endl;
			//}

			std::stringstream residue_key;
			residue_key << pose.pdb_info()->chain(ii) << "," << pose.residue(ii).name3() << "," << pdb_resid;
			// certain structures have extra residues with repeated numbering, so-called insertion codes
			// deal with those cases with the following extra conditional
			if ( icode != ' ' ) {
				residue_key	<< icode;
			}

			float hp_ASA = -123.4;
			if ( res_to_calculated_hp_ASA.find( residue_key.str() ) != res_to_calculated_hp_ASA.end() )
				hp_ASA = res_to_calculated_hp_ASA.find( residue_key.str() )->second;
			else {
				// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
				if ( pose.residue(ii).name3() == "MET" ) {
					residue_key.str("");
					residue_key << pose.pdb_info()->chain(ii) << ",MSE," << pdb_resid;

					if ( res_to_calculated_hp_ASA.find( residue_key.str() ) != res_to_calculated_hp_ASA.end() ) {
						hp_ASA = res_to_calculated_hp_ASA.find( residue_key.str() )->second;
					} else {
						TR << "Residue " << pose.residue( ii ).name3() << " " << ii
							<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
							<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
						continue;
						//exit(1);
					}
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pdb_resid
						<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
						<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
					continue;
					//exit(1);
				}
			}

			// store the hp_ASA into the right vector, i.e. based on the number of neighbors this residue has
			if ( countNeighbors <= 10 ) {
				hp_ASA_by_aa_type_nbs_1_to_10[ pose.residue(ii).aa() ].push_back( hp_ASA );
			} else if ( countNeighbors >= 11 && countNeighbors <= 13 ) {
				hp_ASA_by_aa_type_nbs_11_to_13[ pose.residue(ii).aa() ].push_back( hp_ASA );
			} else if ( countNeighbors >= 14 && countNeighbors <= 16 ) {
				hp_ASA_by_aa_type_nbs_14_to_16[ pose.residue(ii).aa() ].push_back( hp_ASA );
			} else {
				// must have 17, 18, 19, or 20 neighbors
				hp_ASA_by_aa_type_nbs_17_to_20[ pose.residue(ii).aa() ].push_back( hp_ASA );
			}

			#ifdef FILE_DEBUG
				TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << countNeighbors << ", hp_ASA: " << hp_ASA
					<< "; pdb: " << pose.pdb_info()->name() << std::endl;
			#endif

		}

	} // end for loop over multiple input pdb files

	// output the data into separate files for each residue type
	std::stringstream filename;
	std::ofstream hp_ASA_file;

	for ( Size ii=1; ii <= hp_ASA_by_aa_type_nbs_1_to_10.size(); ++ii ) {
		filename << "hp_ASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs1-10.txt";
		if ( trytry_ofstream_open( hp_ASA_file, filename.str(), std::ios::out ) ) {
			hp_ASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs1-10" << std::endl;
			for ( Size jj=1; jj <= hp_ASA_by_aa_type_nbs_1_to_10[ ii ].size(); ++jj ) {
				hp_ASA_file << hp_ASA_by_aa_type_nbs_1_to_10[ ii ][ jj ] << std::endl;
			}
			hp_ASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

	for ( Size ii=1; ii <= hp_ASA_by_aa_type_nbs_11_to_13.size(); ++ii ) {
		filename << "hp_ASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs11-13.txt";
		if ( trytry_ofstream_open( hp_ASA_file, filename.str(), std::ios::out ) ) {
			hp_ASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs11-13" << std::endl;
			for ( Size jj=1; jj <= hp_ASA_by_aa_type_nbs_11_to_13[ ii ].size(); ++jj ) {
				hp_ASA_file << hp_ASA_by_aa_type_nbs_11_to_13[ ii ][ jj ] << std::endl;
			}
			hp_ASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}


	for ( Size ii=1; ii <= hp_ASA_by_aa_type_nbs_14_to_16.size(); ++ii ) {
		filename << "hp_ASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs14-16.txt";
		if ( trytry_ofstream_open( hp_ASA_file, filename.str(), std::ios::out ) ) {
			hp_ASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs14-16" << std::endl;
			for ( Size jj=1; jj <= hp_ASA_by_aa_type_nbs_14_to_16[ ii ].size(); ++jj ) {
				hp_ASA_file << hp_ASA_by_aa_type_nbs_14_to_16[ ii ][ jj ] << std::endl;
			}
			hp_ASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

	for ( Size ii=1; ii <= hp_ASA_by_aa_type_nbs_17_to_20.size(); ++ii ) {
		filename << "hp_ASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs17-20.txt";
		if ( trytry_ofstream_open( hp_ASA_file, filename.str(), std::ios::out ) ) {
			hp_ASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs17-20" << std::endl;
			for ( Size jj=1; jj <= hp_ASA_by_aa_type_nbs_17_to_20[ ii ].size(); ++jj ) {
				hp_ASA_file << hp_ASA_by_aa_type_nbs_17_to_20[ ii ][ jj ] << std::endl;
			}
			hp_ASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

} // end calculate_hydrophobic_accessible_surface_area


/// @brief a quick protocol written to tally the amount of hASA exposed by each residue type in a certain environment

/// @details
// Uses the same reformatted, residue-level accessibility files that were generated by using grep calls on the output files
// generated by NACCESS.  See comments in calculate_hydrophobic_accessible_surface_area for more information.
//
void
calculate_hASA_by_type_and_exposure( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has, broken up by res type and exposure amount
	// we might need to add pseudocounts to this in case certain aa types are not seen at a given exposure level
	// in that case, we could init the vector to 1 in each case
	//utility::vector1< utility::vector1< utility::vector1< core::Real > > > hASA_by_type_and_NV_value( 20, utility::vector1< utility::vector1< Real > >( chemical::num_canonical_aas ) );
	utility::vector1< utility::vector1< utility::vector1< core::Real > > > hASA_by_type_and_nbcount( 20, utility::vector1< utility::vector1< Real > >( chemical::num_canonical_aas ) );

	// iterate through all the structures - do something to them
	for (std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
		// First, don't include the 3-residue chain termini. That could skew the distributions.  People add purification
		// tags and extra residues on proteins all the time.
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		scorefxn->( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hASA;
		read_NACCESS_asa_file( naccess_asa_filename.str(), res_to_calculated_hASA );

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) { continue; }

			/*
			// Next, we have to determine the NV value

			//// the following lifted straight from NVscore.cc written by Sam DeLuca, with extra comments ////
			//rj lbound defaults to 3.3 and ubound defaults to 11.1. see Durham et al, JMolModel, 2009
			Real lBound = 3.3;  // core::options::option[ core::options::OptionKeys::score::NV_lbound]();
			Real uBound = 11.1; // core::options::option[ core::options::OptionKeys::score::NV_ubound]();

			Real neighborCount = 0;
			Vector neighborVectSum( 0, 0, 0 );

			// use the coordinates of residue neighbor atom for all calcuations
			Vector currentVect( pose.residue(ii).nbr_atom_xyz() );

			pose::Pose pose2(pose);
			//rj iterate over all positions to find neighbors
			for ( conformation::ResidueOPs::iterator iter = pose2.res_begin(); iter != pose2.res_end() ; ++iter) {

				// get the residue to compare to the current residue rsd
				conformation::Residue compRes(**iter);
				if ( pose.residue(ii).seqpos() == compRes.seqpos() )
					continue;

				Vector compVect( compRes.nbr_atom_xyz() );
				Vector::Value dist = currentVect.distance(compVect);

				// get the weighted neighbor count
				Real weight = neighborWeight( dist, lBound, uBound );

				//calculate the weighted neighbor vector for this pair and sum
				Vector weightedVector = ( (compVect-currentVect) / dist) * weight;
				neighborCount += weight;
				neighborVectSum += weightedVector;

			}

			Vector avgSum = neighborVectSum / neighborCount;
			// neighbor vector score is the norm of the average sum of all neighbor vectors
			Real NV_value = avgSum.norm();

			//rj this value can be pretty continuous but we want to lump it into bins of size 0.05.
			//rj we can use floor and casting to get it to the right bin size. this wouldn't work for negative numbers
			// NV_value = 0.32878
			// NV_value * 20 = 6.575600000000006
			// NV_value + 0.5 = 7.0756000000006
			// floor( NV_value ) = 7.0
			// NV_value / 20 = 0.35
			NV_value = floor( NV_value * 20 + 0.5 ) / 20;
			// So the smallest value we can achieve is 0.0. Then 0.05. Max is 1.0.  If we multiply by 20, we can get the vector index.
			// Except in the case where we get 1.0 for the NV value.  Then, we need to use index 20, not 21.
			int NV_value_index = int( NV_value * 20 ) + 1;
			if ( NV_value_index > 20 ) { NV_value_index = 20; } // reset to 20 for those rare cases where NV value is 1.0 exactly.

			// got the NV value...

			// do st similar for the neighbor count. it's going to vary from 0 to 20? 30?  let's cut it off at 20, for 20 bins?
			// if we were using negative numbers, this would be better: static_cast<int>(x + x > 0.0 ? +0.5: -0.5);
			int roundedNeighborCount = int( neighborCount + 0.5 );
			if ( roundedNeighborCount < 1 ) { roundedNeighborCount = 1; }
			if ( roundedNeighborCount > 20 ) { roundedNeighborCount = 20; }
			*/

			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			Size neighborCount = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( neighborCount > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// something else to be careful about is that mini likes to renumber the residues that it reads from a PDB file
			// so if the PDB starts at something other than 1, we need to translate the mini number to the PDB number

			// convert the PDB resid to the pose resid in case the PDB numbering is screwy
			int pdb_resid = pose.pdb_info()->number(ii);
			char icode = pose.pdb_info()->icode(ii);

			//if ( pdb_resid != ii ) {
			//	TR << "Converted mini res num " << ii << " to PDB res num " << pdb_resid << std::endl;
			//}

			std::stringstream residue_key;
			residue_key << pose.pdb_info()->chain(ii) << "," << pose.residue(ii).name3() << "," << pdb_resid;

			// certain structures have extra residues with repeated numbering, so-called insertion codes
			// deal with those cases with the following extra conditional
			if ( icode != ' ' ) { residue_key << icode; }

			float hASA = -12345.6;
			if ( res_to_calculated_hASA.find( residue_key.str() ) != res_to_calculated_hASA.end() )
				hASA = res_to_calculated_hASA.find( residue_key.str() )->second;

			else {
				// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
				if ( pose.residue(ii).name3() == "MET" ) {
					residue_key.str("");
					residue_key << pose.pdb_info()->chain(ii) << ",MSE," << pdb_resid;

					if ( res_to_calculated_hASA.find( residue_key.str() ) != res_to_calculated_hASA.end() ) {
						hASA = res_to_calculated_hASA.find( residue_key.str() )->second;
					} else {
						TR << "Residue " << pose.residue( ii ).name3() << " " << ii
							<< " is surface-exposed, with " << neighborCount << " neighbors, but has no hASA info; residue key: "
							<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
						exit(1);
					}
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pdb_resid
						<< ") is surface-exposed, with " << neighborCount << " neighbors, but has no hASA info; residue key: "
						<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
					exit(1);
				}
			}

			//TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << neighborCount
			//	<< ", hASA: " << hASA
			//	<< ", NV_value: " << NV_value << ", NV_value_index: " << NV_value_index
			//	<< "; pdb: " << pose.pdb_info()->name() << std::endl;

			// store the hASA into the right vector, i.e. based on the number of neighbors this residue has
			//hASA_by_type_and_NV_value[ NV_value_index ][ pose.residue(ii).aa() ].push_back( hASA );
			hASA_by_type_and_nbcount[ neighborCount ][ pose.residue(ii).aa() ].push_back( hASA );

		}

	} // end for loop over multiple input pdb files

	// output the data into separate files for each residue type
	std::stringstream filename;
	std::ofstream hASA_file;

	//for ( Size ii=1; ii <= hASA_by_type_and_NV_value.size(); ++ii ) {
	//	for ( Size jj=1; jj <= hASA_by_type_and_NV_value[ ii ].size(); ++jj ) {
	//		filename << "hASA_by_NV" << ii << "." << chemical::name_from_aa( (chemical::AA)jj ) << ".txt";
	//		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
	//			hASA_file << chemical::name_from_aa( (chemical::AA)jj ) << std::endl;
	//			for ( Size kk=1; kk <= (hASA_by_type_and_NV_value[ ii ][ jj ]).size(); ++kk ) {
	//				hASA_file << hASA_by_type_and_NV_value[ ii ][ jj ][ kk ] << std::endl;
	//			}
	//			hASA_file.close();
	//		} else {
	//			TR << "Could not open filename " << filename.str() << std::endl;
	//		}
	//		filename.str("");
	//	}
	//}

	for ( Size ii=1; ii <= hASA_by_type_and_nbcount.size(); ++ii ) {
		for ( Size jj=1; jj <= hASA_by_type_and_nbcount[ ii ].size(); ++jj ) {
			filename << "hASA_by_nbcount" << ii << "." << chemical::name_from_aa( (chemical::AA)jj ) << ".txt";
			if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
				hASA_file << chemical::name_from_aa( (chemical::AA)jj ) << std::endl;
				for ( Size kk=1; kk <= (hASA_by_type_and_nbcount[ ii ][ jj ]).size(); ++kk ) {
					hASA_file << hASA_by_type_and_nbcount[ ii ][ jj ][ kk ] << std::endl;
				}
				hASA_file.close();
			} else {
				TR << "Could not open filename " << filename.str() << std::endl;
			}
			filename.str("");
		}
	}


} // end calculate_hASA_by_type_and_exposure


void
calculate_hASA_by_type_and_attractiveE( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has, broken up by res type and exposure amount
	// we might need to add pseudocounts to this in case certain aa types are not seen at a given exposure level
	// in that case, we could init the vector to 1 in each case
	//utility::vector1< utility::vector1< utility::vector1< core::Real > > > hASA_by_type_and_atrE( 32, utility::vector1< utility::vector1< Real > >( chemical::num_canonical_aas ) );
	utility::vector1< utility::vector1< utility::vector1< core::Real > > > hASA_by_type_and_atrE( 16, utility::vector1< utility::vector1< Real > >( chemical::num_canonical_aas ) );

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		scorefxn->score( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hASA;
		read_NACCESS_asa_file( naccess_asa_filename.str(), res_to_calculated_hASA );

		// First, don't include the 3-residue chain termini. That could skew the distributions.  People add purification
		// tags and extra residues on proteins all the time.
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) { continue; }

			// Next, we have to determine the NV value

			////// the following lifted straight from NVscore.cc written by Sam DeLuca, with extra comments /////
			//rj lbound defaults to 3.3 and ubound defaults to 11.1. see Durham et al, JMolModel, 2009
			Real lBound = 3.3;  // core::options::option[ core::options::OptionKeys::score::NV_lbound]();
			Real uBound = 11.1; // core::options::option[ core::options::OptionKeys::score::NV_ubound]();

			Real neighborCount = 0;
			Vector neighborVectSum( 0, 0, 0 );

			// use the coordinates of residue neighbor atom for all calcuations
			Vector currentVect( pose.residue(ii).nbr_atom_xyz() );

			pose::Pose pose2(pose);
			//rj iterate over all positions to find neighbors
			//for ( conformation::ResidueOPs::iterator iter = pose2.res_begin(); iter != pose2.res_end() ; ++iter) {
			for ( Size ii_res = 1; ii_res <= pose2.size(); ++ii_res ) {

				// get the residue to compare to the current residue rsd
				//conformation::Residue compRes(**iter);
				conformation::Residue const & compRes = pose2.residue( ii_res );
				if ( pose.residue(ii).seqpos() == compRes.seqpos() )
					continue;

				Vector compVect( compRes.nbr_atom_xyz() );
				Vector::Value dist = currentVect.distance(compVect);

				// get the weighted neighbor count
				Real weight = neighborWeight( dist, lBound, uBound );

				//calculate the weighted neighbor vector for this pair and sum
				Vector weightedVector = ( (compVect-currentVect) / dist) * weight;
				neighborCount += weight;
				neighborVectSum += weightedVector;

			}

			Vector avgSum = neighborVectSum / neighborCount;
			// neighbor vector score is the norm of the average sum of all neighbor vectors
			Real NV_value = avgSum.norm();

			//rj this value can be pretty continuous but we want to lump it into bins of size 0.05.
			//rj we can use floor and casting to get it to the right bin size. this wouldn't work for negative numbers
			// NV_value = 0.32878
			// NV_value * 20 = 6.575600000000006
			// NV_value + 0.5 = 7.0756000000006
			// floor( NV_value ) = 7.0
			// NV_value / 20 = 0.35
			NV_value = floor( NV_value * 20 + 0.5 ) / 20;
			// So the smallest value we can achieve is 0.0. Then 0.05. Max is 1.0.  If we multiply by 20, we can get the vector index.
			// Except in the case where we get 1.0 for the NV value.  Then, we need to use index 20, not 21.
			int NV_value_index = int( NV_value * 20 ) + 1;
			if ( NV_value_index > 20 ) { NV_value_index = 20; } // reset to 20 for those rare cases where NV value is 1.0 exactly.

			// got the NV value...

			// do st similar for the neighbor count. it's going to vary from 0 to 20? 30?  let's cut it off at 17, for 16 bins.
			// if we were using negative numbers, this would be better: static_cast<int>(x + x > 0.0 ? +0.5: -0.5);
			int roundedNeighborCount = int( neighborCount + 0.5 );
			if ( roundedNeighborCount < 1 ) { roundedNeighborCount = 1; }
			if ( roundedNeighborCount > 16 ) { roundedNeighborCount = 16; }


			// round this residues fa_atr energy to the nearest 0.5 kcal
			// off a quick test, fa_atr residue energies range from 0 to -15. So let's use a size of 32.
			float atrE = (pose.energies().residue_total_energies( ii ))[ core::scoring::fa_atr ];

			//float rounded_atrE = ceil( atrE * 2 - 0.5 ) / 2;
			//-2.2498, -2.0
			// -2.495, -2.5
			// -2.510, -2.5,
			// -2.752, -3.0
			// -2.7499, -2.5
			//int atrE_index = int( rounded_atrE * -2 );
			//if ( atrE_index < 1 ) { atrE_index = 1; }
			//if ( atrE_index > 32 ) { atrE_index = 32; }
			// -0.346, -0.5, 1
			// -0.2499, 0, 0
			// -13.5055, -13.5, 27
			// -13.91, -14, 28

			float rounded_atrE = ceil( atrE - 0.5 );
			// -2.2498, -2.0
			// -2.495, -2.0
			// -2.510, -3.0,
			// -2.752, -3.0
			// -2.7499, -3.0
			int atrE_index = int( rounded_atrE * -1 );
			if ( atrE_index < 1 ) { atrE_index = 1; }
			if ( atrE_index > 16 ) { atrE_index = 16; }

			// something else to be careful about is that mini likes to renumber the residues that it reads from a PDB file
			// so if the PDB starts at something other than 1, we need to translate the mini number to the PDB number

			// convert the PDB resid to the pose resid in case the PDB numbering is screwy
			int pdb_resid = pose.pdb_info()->number(ii);
			char icode = pose.pdb_info()->icode(ii);

			std::stringstream residue_key;
			residue_key << pose.pdb_info()->chain(ii) << "," << pose.residue(ii).name3() << "," << pdb_resid;

			// certain structures have extra residues with repeated numbering, so-called insertion codes
			// deal with those cases with the following extra conditional
			if ( icode != ' ' ) { residue_key << icode; }

			float hASA = -12345.6;
			if ( res_to_calculated_hASA.find( residue_key.str() ) != res_to_calculated_hASA.end() )
				hASA = res_to_calculated_hASA.find( residue_key.str() )->second;

			else {
				// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
				if ( pose.residue(ii).name3() == "MET" ) {
					residue_key.str("");
					residue_key << pose.pdb_info()->chain(ii) << ",MSE," << pdb_resid;

					if ( res_to_calculated_hASA.find( residue_key.str() ) != res_to_calculated_hASA.end() ) {
						hASA = res_to_calculated_hASA.find( residue_key.str() )->second;
					} else {
						TR << "Residue " << pose.residue( ii ).name3() << " " << ii
							<< " is surface-exposed, with " << neighborCount << " neighbors, but has no hASA info; residue key: "
							<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
						exit(1);
					}
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pdb_resid
						<< ") is surface-exposed, with " << neighborCount << " neighbors, but has no hASA info; residue key: "
						<< residue_key.str() << ", pdb: " << pose.pdb_info()->name() << std::endl;
					exit(1);
				}
			}

			//TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << neighborCount
			//	<< ", roundedNeighborCount: " << roundedNeighborCount << ", hASA: " << hASA
			//	<< ", NV_value: " << NV_value << ", NV_value_index: " << NV_value_index
			//	<< ", atrE: " << atrE << ", rounded_atrE: " << rounded_atrE << ", atrE_index: " << atrE_index
			//	<< "; pdb: " << pose.pdb_info()->name() << std::endl;

			// store the hASA into the right vector, i.e. based on the fa_atr energy this residue has
			hASA_by_type_and_atrE[ atrE_index ][ pose.residue(ii).aa() ].push_back( hASA );

		}

	} // end for loop over multiple input pdb files

	// output the data into separate files for each residue type
	std::stringstream filename;
	std::ofstream hASA_file;

	for ( Size ii=1; ii <= hASA_by_type_and_atrE.size(); ++ii ) {
		for ( Size jj=1; jj <= hASA_by_type_and_atrE[ ii ].size(); ++jj ) {
			filename << "hASA_by_atrE" << ii << "." << chemical::name_from_aa( (chemical::AA)jj ) << ".txt";
			if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
				hASA_file << chemical::name_from_aa( (chemical::AA)jj ) << std::endl;
				for ( Size kk=1; kk <= (hASA_by_type_and_atrE[ ii ][ jj ]).size(); ++kk ) {
					hASA_file << hASA_by_type_and_atrE[ ii ][ jj ][ kk ] << std::endl;
				}
				hASA_file.close();
			} else {
				TR << "Could not open filename " << filename.str() << std::endl;
			}
			filename.str("");
		}
	}


} // end calculate_hASA_by_type_and_attractiveE


/// @brief helper function for function below
std::string
get_residue_key( core::pose::Pose & pose, core::Size pose_resid, bool mse ) {

	// something to be careful about is that mini likes to renumber the residues that it reads from a PDB file
	// so if the PDB starts at something other than 1, we need to translate the mini number to the PDB number

	// convert the PDB resid to the pose resid in case the PDB numbering is screwy
	int pdb_resid = pose.pdb_info()->number( pose_resid );
	char icode = pose.pdb_info()->icode( pose_resid );

	//if ( pdb_resid != pose_resid ) {
	//	TR << "Converted mini res num " << pose_resid << " to PDB res num " << pdb_resid << std::endl;
	//}

	std::stringstream residue_key;
	residue_key << pose.pdb_info()->chain(pose_resid) << ",";
	if ( mse ) {
		residue_key << "MSE,";
	} else {
		residue_key << pose.residue(pose_resid).name3() << ",";
	}
	residue_key << pdb_resid;

	// certain structures have extra residues with repeated numbering, so-called insertion codes
	// deal with those cases with the following extra conditional
	if ( icode != ' ' ) {
		residue_key	<< icode;
	}

	return residue_key.str();
}


/// @details
/// Same as function above but instead of using the average amount of hASA a residue exposes in a certain environment, go back
/// and read the NACCESS files to see how much surface area exactly the neighbors and that residue expose. Doing this to get an
/// idea of the amount of error we introduce by using average hASA instead of exact hASA.
/// After going through all PDBs and all residues, it will print out the hASA in that sphere and I will create a distribution of this in JMP.
///
void
calculate_total_hASA_within_distance_exact_hASA_values( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< Real > hASA_within_10A;


	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		scorefxn->score( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hASA;
		read_NACCESS_asa_file( naccess_asa_filename.str(), res_to_calculated_hASA );

		// First, we have to figure out which residues are surface-exposed
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		Real patchArea = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hpASA (based on its
		// neighbor count) and all neighboring surface-exposed residues hpASA's to the total.  Save that final sum in the
		// hASA_within_10A vector.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) {
				continue;
			}

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// passed the surface-exposed check...

			// this converts a pose residue index into the PDB residue number (in case the PDB doesn't start at 1)
			std::string key = get_residue_key( pose, ii, false /*not a MET*/ );

			// reset the area size for every surface-exposed residue
			patchArea = 0.0;

			// don't add in the hASA contribution if the residue is polar
			if ( ! pose.residue(ii).is_polar() ) {

				if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
					patchArea += res_to_calculated_hASA.find( key )->second;
				else {
					// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
					if ( pose.residue( ii ).name3() == "MET" ) {
						key = get_residue_key( pose, ii, true /*useMSE*/ );
						if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
							patchArea += res_to_calculated_hASA.find( key )->second;
						} else {
							TR << "Residue " << pose.residue( ii ).name3() << " " << ii
								<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
								<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
							exit(1);
						}
					} else {
						TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pose.pdb_info()->number( ii )
							<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
							<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
						exit(1);
					}
				}
			}

			// add in every neighbors hASA too
			// for every Edge in the neighbor graph, figure out if that residue is surface exposed

			//TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " ( " << patchArea << " ) include " << std::endl;
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( !(pose.residue(jj).is_protein()) ) {
					continue;
				}

				if ( pose.residue(jj).is_polar() ) {
					continue;
				}

				// BUG FIX: Don't use the same residue key that we got before. Otherwise we just keep adding residue ii's hASA to the total.
				// Instead, make a new key for each neighboring residue, and add that residue's hASA to the total.
				key = get_residue_key( pose, jj, false /*not MET*/ );

				//TR << pose.residue(jj).name3() << " " << jj << ": ";

				Real area = 0.0;
				countNeighbors = tenA_neighbor_graph.get_node(jj)->num_neighbors_counting_self();
				if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
					//TR << "Not SE" << std::endl;
					continue;
				} else {
					if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
						area = res_to_calculated_hASA.find( key )->second;
					else {
						// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
						if ( pose.residue( jj ).name3() == "MET" ) {
							key = get_residue_key( pose, jj, true /*useMSE*/ );
							if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
								area = res_to_calculated_hASA.find( key )->second;
							} else {
								TR << "Residue " << pose.residue( jj ).name3() << " " << jj
									<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
									<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
								exit(1);
							}
						} else {
							TR << "Residue " << pose.residue( jj ).name3() << " " << jj << " (PDB: " << pose.pdb_info()->number( jj )
								<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
								<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
							exit(1);
						}
					}
					patchArea += area;
				}
				//TR << area << std::endl;
			}
			//TR << "\n" << std::endl;

			// now that we know how much exposed hASA this residue has, save it somewhere
			if ( patchArea != 0.0 ) {
				hASA_within_10A.push_back( patchArea );

				//if ( patchArea < 1.0 ) {
				//	TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pose.pdb_info()->number( ii )
				//		<< ") has very small patchArea. residue key: " << key << ", pdb: " << pose.pdb_info()->name();
				//	TR << ", nbs: ";
				//	for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				//		eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {
				//		Size jj = (*eli)->get_other_ind(ii);
				//		TR << pose.residue( jj ).name3() << "-" << jj << "(" << pose.pdb_info()->number( jj ) << "), ";
				//	}
				//	TR << std::endl;
				//}
			}

		} // end loop over all residues

	} // end for loop over multiple input pdb files


	// output the data into separate files for each residue type
	std::ofstream patches_file;
	if ( trytry_ofstream_open( patches_file, "hASA_patch_sizes.using_NACCESS_values.txt", std::ios::out ) ) {
		for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
			patches_file << hASA_within_10A[ii] << std::endl;
		}
		patches_file.close();
	} else {
		TR << "Could not open filename hASA_patch_sizes.txt" << std::endl;
	}


} // end calculate_total_hASA_within_distance_exact_hASA_values


/// @details
/// Same as function above but instead of using just hydrophobic neighbors, uses all surface-exposed neighbors.
///
void
calculate_total_hASA_within_distance_exact_hASA_values_allnbs( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< Real > hASA_within_10A;


	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		scorefxn->( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hASA;
		read_NACCESS_asa_file( naccess_asa_filename.str(), res_to_calculated_hASA );

		// First, we have to figure out which residues are surface-exposed
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		Real patchArea = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hpASA (based on its
		// neighbor count) and all neighboring surface-exposed residues hpASA's to the total.  Save that final sum in the
		// hASA_within_10A vector.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) {
				continue;
			}

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// passed the surface-exposed check...

			// this converts a pose residue index into the PDB residue number (in case the PDB doesn't start at 1)
			std::string key = get_residue_key( pose, ii, false /*not a MET*/ );

			// reset the area size for every surface-exposed residue
			patchArea = 0.0;

			// add in the hASA contribution regardless of whether the residue is polar or hydrophobic
			if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
				patchArea += res_to_calculated_hASA.find( key )->second;
			else {
				// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
				if ( pose.residue( ii ).name3() == "MET" ) {
					key = get_residue_key( pose, ii, true /*useMSE*/ );
					if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
						patchArea += res_to_calculated_hASA.find( key )->second;
					} else {
						TR << "Residue " << pose.residue( ii ).name3() << " " << ii
							<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
							<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
						exit(1);
					}
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pose.pdb_info()->number( ii )
						<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
						<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
					exit(1);
				}
			}

			// add in every neighbors hASA too
			// for every Edge in the neighbor graph, figure out if that residue is surface exposed

			TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " ( " << patchArea << " ) include " << std::endl;
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( !(pose.residue(jj).is_protein()) ) {
					continue;
				}

				// BUG FIX: Don't use the same residue key that we got before. Otherwise we just keep adding residue ii's hASA to the total.
				// Instead, make a new key for each neighboring residue, and add that residue's hASA to the total.
				key = get_residue_key( pose, jj, false /*not MET*/ );

				TR << pose.residue(jj).name3() << " " << jj << ": ";

				Real area = 0.0;
				countNeighbors = tenA_neighbor_graph.get_node(jj)->num_neighbors_counting_self();
				if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
					TR << "Not SE" << std::endl;
					continue;
				} else {
					if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
						area = res_to_calculated_hASA.find( key )->second;
					else {
						// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
						if ( pose.residue( jj ).name3() == "MET" ) {
							key = get_residue_key( pose, jj, true /*useMSE*/ );
							if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
								area = res_to_calculated_hASA.find( key )->second;
							} else {
								TR << "Residue " << pose.residue( jj ).name3() << " " << jj
									<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
									<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
								exit(1);
							}
						} else {
							TR << "Residue " << pose.residue( jj ).name3() << " " << jj << " (PDB: " << pose.pdb_info()->number( jj )
								<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
								<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
							exit(1);
						}
					}
					patchArea += area;
				}
				TR << area << std::endl;
			}
			TR << "\n" << std::endl;

			// now that we know how much exposed hASA this residue has, save it somewhere
			if ( patchArea != 0.0 ) {
				hASA_within_10A.push_back( patchArea );

				//if ( patchArea < 1.0 ) {
				//	TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pose.pdb_info()->number( ii )
				//		<< ") has very small patchArea. residue key: " << key << ", pdb: " << pose.pdb_info()->name();
				//	TR << ", nbs: ";
				//	for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				//		eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {
				//		Size jj = (*eli)->get_other_ind(ii);
				//		TR << pose.residue( jj ).name3() << "-" << jj << "(" << pose.pdb_info()->number( jj ) << "), ";
				//	}
				//	TR << std::endl;
				//}
			}

		} // end loop over all residues

	} // end for loop over multiple input pdb files


	// output the data into separate files for each residue type
	std::ofstream patches_file;
	if ( trytry_ofstream_open( patches_file, "hASA_patch_sizes.using_NACCESS_values.allnbs.txt", std::ios::out ) ) {
		for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
			patches_file << hASA_within_10A[ii] << std::endl;
		}
		patches_file.close();
	} else {
		TR << "Could not open filename hASA_patch_sizes.using_NACCESS_values.allnbs.txt" << std::endl;
	}


} // end calculate_total_hASA_within_distance_exact_hASA_values_allnbs


/// @details
/// Same as function above but instead of using just hydrophobic neighbors, uses all neighbors, surface exposed or not.
/// This function makes multiple score distributions, based on the number of neighbors a position has.
///
void
calculate_total_hASA_within_distance_exact_hASA_values_allnbs_exposedornot_conditionalonnumnbs( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< std::pair< Size, Real > > hASA_within_10A;

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		// check to make sure NACCESS data for this file exists, too
		std::stringstream naccess_asa_filename;
		naccess_asa_filename << utility::file::file_basename( *pdb ) << ".rsa.reformatted";
		if ( !utility::file::file_exists( naccess_asa_filename.str() ) ) {
			std::cerr << "NACCESS asa file for pdb " << *pdb << ", '" << naccess_asa_filename.str() << "' not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		scorefxn->( pose );

		// may as well open up the asa file and read it into memory so we don't have to open and reopen for each residue
		std::map< std::string, Real > res_to_calculated_hASA;
		read_NACCESS_asa_file( naccess_asa_filename.str(), res_to_calculated_hASA );

		Real patchArea = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hASA (based on its
		// neighbor count) and all neighboring surface-exposed residues hASA's to the total.  Save that final sum in the
		// hASA_within_10A vector.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() )
				continue;

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			// only calculate the hpatch data for residues with fewer than 20nbs. but do include all neighboring residues
			// not just the surface exposed ones.
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF )
				continue;

			// passed the surface-exposed check...

			// this converts a pose residue index into the PDB residue number (in case the PDB doesn't start at 1)
			std::string key = get_residue_key( pose, ii, false /*not a MET*/ );

			// reset the area size for every surface-exposed residue
			patchArea = 0.0;

			// add in the hASA contribution regardless of whether the residue is polar or hydrophobic
			if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
				patchArea += res_to_calculated_hASA.find( key )->second;
			else {
				// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
				if ( pose.residue( ii ).name3() == "MET" ) {
					key = get_residue_key( pose, ii, true /*useMSE*/ );
					if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
						patchArea += res_to_calculated_hASA.find( key )->second;
					} else {
						//TR << "Residue " << pose.residue( ii ).name3() << " " << ii
						//	<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
						//	<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
						//exit(1);
						continue;
					}
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << " (PDB: " << pose.pdb_info()->number( ii )
						<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
						<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
					//exit(1);
				}
			}

			// add in every neighbors hASA too

			//TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " ( " << patchArea << " ) include " << std::endl;
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( !(pose.residue(jj).is_protein()) ) { continue; }

				// BUG FIX: Don't use the same residue key that we got before. Otherwise we just keep adding residue ii's hASA to the total.
				// Instead, make a new key for each neighboring residue, and add that residue's hASA to the total.
				key = get_residue_key( pose, jj, false /*not MET*/ );
				//TR << pose.residue(jj).name3() << " " << jj << ": ";

				Real area = 0.0;
				if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() )
					area = res_to_calculated_hASA.find( key )->second;
				else {
					// maybe the residue is a seleno-MET; look for that key and if that's not found then we give up
					if ( pose.residue( jj ).name3() == "MET" ) {
						key = get_residue_key( pose, jj, true /*useMSE*/ );
						if ( res_to_calculated_hASA.find( key ) != res_to_calculated_hASA.end() ) {
							area = res_to_calculated_hASA.find( key )->second;
						} else {
							//TR << "Residue " << pose.residue( jj ).name3() << " " << jj
							//	<< " is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
							//	<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
							//exit(1);
							continue;
						}
					} else {
						TR << "Residue " << pose.residue( jj ).name3() << " " << jj << " (PDB: " << pose.pdb_info()->number( jj )
							<< ") is surface-exposed, with " << countNeighbors << " neighbors, but has no hASA info; residue key: "
							<< key << ", pdb: " << pose.pdb_info()->name() << std::endl;
						//exit(1);
					}
				}
				patchArea += area;
			}

			// now that we know how much exposed hASA this residue has, save it somewhere
			hASA_within_10A.push_back( std::make_pair( countNeighbors, patchArea ) );

		} // end loop over all residues

	} // end for loop over multiple input pdb files


	// output the data into separate files for each residue type
	std::ofstream patches_file;
	if ( trytry_ofstream_open( patches_file, "hASA_patch_sizes.using_NACCESS_values.allnbs.txt", std::ios::out ) ) {
		for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
			patches_file << hASA_within_10A[ii].first << '\t' << hASA_within_10A[ii].second << std::endl;
		}
		patches_file.close();
	} else {
		TR << "Could not open filename hASA_patch_sizes.using_NACCESS_values.allnbs.txt" << std::endl;
	}


} // end calculate_total_hASA_within_distance_exact_hASA_values_allnbs_exposedornot_conditionalonnumnbs


/// @details
/// Same as function above but instead of using just hydrophobic neighbors, uses all neighbors, surface exposed or not.
/// This function makes multiple score distributions, based on the number of neighbors a position has.
/// Same as above but instead of using NACCESS, uses the SurfacePotential to get hASA.
///
void
calculate_total_hASA_within_distance_miniSASAvalues_allnbs_exposedornot_conditionalonnumnbs( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< utility::vector1< Real > > hASA_within_10A( 20 );

	// iterate through all the structures - do something to them
	for ( std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb ) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		scorefxn->score( pose );

		utility::vector1< core::Real > residue_sasa( pose.size(), 0.0 );
		utility::vector1< core::Real > residue_hsasa( pose.size(), 0.0 ); // hydrophobic SASA only
		core::Real total_hydrophobic_sasa = core::scoring::calc_per_res_hydrophobic_sasa( pose, residue_sasa, residue_hsasa, 1.4 /* probe radius */ );
		#ifdef FILE_DEBUG
			TR << "total_hydrophobic_sasa: " << total_hydrophobic_sasa << std::endl;
		#endif

		Real total_area = 0.0;
		Real total_hydrophobic_area = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hASA (based on its
		// neighbor count) and all neighboring surface-exposed residues hASA's to the total.  Save that final sum in the
		// hASA_within_10A vector.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() )
				continue;

			// our definition of surface residue is that the residue has fewer than 20 neighbors
			// only calculate the hpatch data for residues with fewer than 20nbs. but do include all neighboring residues
			// not just the surface exposed ones.
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF )
				continue;

			// passed the surface-exposed check...
			// reset the area size for every surface-exposed residue
			total_area = 0.0;
			total_hydrophobic_area = 0.0;

			// add in the hASA contribution of this residue regardless of whether the residue is polar or hydrophobic
			total_area += residue_sasa[ ii ];
			total_hydrophobic_area += residue_hsasa[ ii ];

			// add in every neighbors hASA too

			//TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " ( " << residue_hsasa[ii] << " ) include " << std::endl;
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( !(pose.residue(jj).is_protein()) )
					continue;

				total_area += residue_sasa[ jj ];
				total_hydrophobic_area += residue_hsasa[ jj ];
			}

			//utility::vector1< Real > values;
			//values.push_back( total_hydrophobic_area );
			//values.push_back( total_area );
			//values.push_back( total_hydrophobic_area / total_area );

			// now that we know how much exposed hASA this residue has, save it somewhere
			//hASA_within_10A[ countNeighbors ].push_back( values );
			hASA_within_10A[ countNeighbors ].push_back( total_hydrophobic_area );

		} // end loop over all residues

	} // end for loop over multiple input pdb files

	// output the data into separate files for each residue type
	std::stringstream filename;

	std::ofstream patches_file;
	for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
		filename << "hASA_patch_sizes.miniSASAvalues.allnbs.txt.nb" << I(2,ii);
		if ( trytry_ofstream_open( patches_file, filename.str(), std::ios::out ) ) {
			patches_file << "nb" << ii << std::endl;
			for ( Size jj=1; jj <= hASA_within_10A[ ii ].size(); ++jj ) {
				patches_file << hASA_within_10A[ ii ][ jj ] << std::endl;
			}
			patches_file.close();
			filename.str("");
		} else {
			TR << "Could not open filename " << filename << std::endl;
		}
	}

} // end calculate_total_hASA_within_distance_miniSASAvalues_allnbs_exposedornot_conditionalonnumnbs


/// @brief helper function for function below
std::string
get_map_key( std::string resname, core::Size count_nbs ) {

	std::stringstream residue_key;
	residue_key << resname;
	if ( count_nbs <= 10 ) {
		residue_key << "-lte10";
	} else if ( count_nbs >= 11 && count_nbs <= 13 ) {
		residue_key << "-lte13";
	} else if ( count_nbs >= 14 && count_nbs <= 16 ) {
		residue_key << "-lte16";
	} else {
		// must have 17, 18, 19 or 20 neighbors
		residue_key << "-lte20";
	}
	return residue_key.str();
}

/// @brief a quick protocol written to calculate hydrophobic ASA on a set of proteins

/// @details
/// The function before this one figures out what the hASA is for each residue type on the surface. This functions uses the values
/// determined in that function. For every surface exposed residue, we look at all the other surface-exposed residues within
/// some distance (e.g. 10A to start with) and how many neighbors THOSE residues have and determine what the total hASA for the
/// surface-exposed residue we're currently on.  After going through all PDBs and all residues, it will print out the hASA
/// in that sphere and I will create a distribution of this in JMP.
///
void
calculate_total_hASA_within_distance_avg_values( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< Real > hASA_within_10A;

	// read in the per-residue hASA means
	// they'll be accessible via this map like: map[ ALA-lte10 ], map[ GLN-lte13 ], map[ TYR-lte16 ]
	std::map< std::string, Real > res_to_average_hASA;

	utility::io::izstream residue_hASA_ifstream;
	basic::database::open( residue_hASA_ifstream, "scoring/score_functions/SurfacePotential/average_hASA_by_res_and_neighbor.txt" );

	std::string restype;
	Real lte10_asa = 0.0;
	Real lte13_asa = 0.0;
	Real lte16_asa = 0.0;

	while ( !residue_hASA_ifstream.eof() ) {
		residue_hASA_ifstream >> restype >> lte10_asa >> lte13_asa >> lte16_asa;

		res_to_average_hASA[ restype + "-lte10" ] = lte10_asa;
		res_to_average_hASA[ restype + "-lte13" ] = lte13_asa;
		res_to_average_hASA[ restype + "-lte16" ] = lte16_asa;
	}

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);

		scorefxn->score( pose );

		// First, we have to figure out which residues are surface-exposed
		// Except, don't include the 3-residue chain termini.  Lots of protein expression systems use
		// MET on the ends which end up on the termini in the crystal structures
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		Real patchArea = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hASA (based on its
		// neighbor count) and all neighboring surface-exposed residues hpASA's to the total.  Save that final sum in the
		// hASA_within_10A vector.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) {
				continue;
			}

			// our definition of surface residue is that the residue has fewer than 16 neighbors
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
				// not a surface exposed residue
				continue;
			}

			// reset the area size for every surface-exposed residue
			patchArea = 0.0;

			// don't add in the hASA contribution if the residue is polar
			std::string key = get_map_key( pose.residue(ii).name3(), countNeighbors );
			if ( ! pose.residue(ii).is_polar() ) {

				if ( res_to_average_hASA.find( key ) != res_to_average_hASA.end() ) {
					patchArea += res_to_average_hASA.find( key )->second;
					//TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << countNeighbors << ", hp_ASA: " << res_to_average_hASA.find( key )->second
					//	<< "; pdb: " << pose.pdb_info()->name() << std::endl;
				} else {
					TR << "Residue " << pose.residue( ii ).name3() << " " << ii << ", nbs: "
						<< countNeighbors << ", no hASA info; residue key: " << key << ", pdb: "
						<< pose.pdb_info()->name() << std::endl;
					exit(1);
				}
			}

			// add in every neighbors hASA too
			// for every Edge in the neighbor graph, figure out if that residue is surface exposed

			//TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " include " << std::endl;
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( ! pose.residue(jj).is_protein() ) { continue; }
				if ( pose.residue(jj).is_polar() ) { continue; }

				//TR << pose.residue(jj).name3() << " " << jj;

				countNeighbors = tenA_neighbor_graph.get_node(jj)->num_neighbors_counting_self();
				if ( countNeighbors > SURFACE_EXPOSED_CUTOFF ) {
					//TR << std::endl;
					continue;
				} else {
					key = get_map_key( pose.residue(jj).name3(), countNeighbors );

					if ( res_to_average_hASA.find( key ) != res_to_average_hASA.end() ) {
						//TR << ": nbs: " << countNeighbors << ", hASA: " << res_to_average_hASA.find( key )->second
						//<< "; pdb: " << pose.pdb_info()->name();
						patchArea += res_to_average_hASA.find( key )->second;
					} else {
						TR << pose.residue(jj).name3() << " " << jj << ", nbs: " << countNeighbors << ", no hASA info; residue key: " << key << std::endl;
						exit(1);
					}
				}
				//TR << std::endl;
			}
			//TR << "\n" << std::endl;

			// now that we know how much exposed hASA this residue has, save it somewhere
			if ( patchArea != 0.0 ) {
				hASA_within_10A.push_back( patchArea );
			}

		} // end loop over all residues

	} // end for loop over multiple input pdb files


	// output the data into separate files for each residue type
	std::ofstream patches_file;
	if ( trytry_ofstream_open( patches_file, "hASA_patch_sizes.txt", std::ios::out ) ) {
		for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
			patches_file << hASA_within_10A[ii] << std::endl;
		}
		patches_file.close();
	} else {
		TR << "Could not open filename hASA_patch_sizes.txt" << std::endl;
	}


} // end calculate_total_hASA_within_distance_avg_values


/// @details
/// Function which determines the amount of hASA a residue exposes given some number of neighbors.
/// Previous versions of this function have used the NACCESS generated output files to determine the hASA. This new version
/// of the function uses the mini SASA calculation machinery to do the same thing.  It will be interesting to see if the
/// two match up.
///
void calculate_hASA_by_type_and_nbcount( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has, broken up by res type
	utility::vector1< utility::vector1< Real > > hASA_by_aa_type_nbs_1_to_10( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hASA_by_aa_type_nbs_11_to_13( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hASA_by_aa_type_nbs_14_to_16( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hASA_by_aa_type_nbs_17_to_20( chemical::num_canonical_aas );
	utility::vector1< utility::vector1< Real > > hASA_by_aa_type_nbs_21_to_24( chemical::num_canonical_aas );

	// iterate through all the structures - do something to them
	for(std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		scorefxn->score( pose );

		utility::vector1< core::Real > residue_sasa( pose.size(), 0.0 );
		utility::vector1< core::Real > residue_hsasa( pose.size(), 0.0 ); // hydrophobic SASA only
		core::Real total_hydrophobic_sasa = core::scoring::calc_per_res_hydrophobic_sasa( pose, residue_sasa, residue_hsasa, 1.4 /* probe radius */, true /* use naccess sasa radii */ );
		#ifdef FILE_DEBUG
			TR << "total_hydrophobic_sasa: " << total_hydrophobic_sasa << std::endl;
		#endif

		// Now iterate over all residues (well, exclude 3 from each end) and figure out which residues are surface-exposed
		float hASA = 0.0;
		for ( Size ii = 4; ii <= pose.size() - 3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() ) continue;

			// our definition of surface residue is that the residue has fewer than 16 or 20 neighbors
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			Size countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
				// don't bother getting stats for such "buried" residues
				continue;
			}

			// passed the surface-exposed check...
			hASA = residue_hsasa[ ii ];

			// store the hASA into the right vector, i.e. based on the number of neighbors this residue has
			if ( countNeighbors <= 10 ) {
				hASA_by_aa_type_nbs_1_to_10[ pose.residue(ii).aa() ].push_back( hASA );
			} else if ( countNeighbors >= 11 && countNeighbors <= 13 ) {
				hASA_by_aa_type_nbs_11_to_13[ pose.residue(ii).aa() ].push_back( hASA );
			} else if ( countNeighbors >= 14 && countNeighbors <= 16 ) {
				hASA_by_aa_type_nbs_14_to_16[ pose.residue(ii).aa() ].push_back( hASA );
			} else if ( countNeighbors >= 17 && countNeighbors <= 20 ) {
				hASA_by_aa_type_nbs_17_to_20[ pose.residue(ii).aa() ].push_back( hASA );
			} else {
				// must have 21, 22, 23, or 24 neighbors
				hASA_by_aa_type_nbs_21_to_24[ pose.residue(ii).aa() ].push_back( hASA );
			}

			#ifdef FILE_DEBUG
				TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << countNeighbors << ", hASA: " << hASA << std::endl;
			#endif

		}

	} // end for loop over multiple input pdb files

	// output the data into separate files for each residue type
	std::stringstream filename;
	std::ofstream hASA_file;

	for ( Size ii=1; ii <= hASA_by_aa_type_nbs_1_to_10.size(); ++ii ) {
		filename << "hASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs1-10.txt";
		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
			hASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs1-10" << std::endl;
			for ( Size jj=1; jj <= hASA_by_aa_type_nbs_1_to_10[ ii ].size(); ++jj ) {
				hASA_file << hASA_by_aa_type_nbs_1_to_10[ ii ][ jj ] << std::endl;
			}
			hASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

	for ( Size ii=1; ii <= hASA_by_aa_type_nbs_11_to_13.size(); ++ii ) {
		filename << "hASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs11-13.txt";
		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
			hASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs11-13" << std::endl;
			for ( Size jj=1; jj <= hASA_by_aa_type_nbs_11_to_13[ ii ].size(); ++jj ) {
				hASA_file << hASA_by_aa_type_nbs_11_to_13[ ii ][ jj ] << std::endl;
			}
			hASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}


	for ( Size ii=1; ii <= hASA_by_aa_type_nbs_14_to_16.size(); ++ii ) {
		filename << "hASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs14-16.txt";
		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
			hASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs14-16" << std::endl;
			for ( Size jj=1; jj <= hASA_by_aa_type_nbs_14_to_16[ ii ].size(); ++jj ) {
				hASA_file << hASA_by_aa_type_nbs_14_to_16[ ii ][ jj ] << std::endl;
			}
			hASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

	for ( Size ii=1; ii <= hASA_by_aa_type_nbs_17_to_20.size(); ++ii ) {
		filename << "hASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs17-20.txt";
		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
			hASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs17-20" << std::endl;
			for ( Size jj=1; jj <= hASA_by_aa_type_nbs_17_to_20[ ii ].size(); ++jj ) {
				hASA_file << hASA_by_aa_type_nbs_17_to_20[ ii ][ jj ] << std::endl;
			}
			hASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

	for ( Size ii=1; ii <= hASA_by_aa_type_nbs_21_to_24.size(); ++ii ) {
		filename << "hASA_by_res." << chemical::name_from_aa( (chemical::AA)ii ) << ".nbs21-24.txt";
		if ( trytry_ofstream_open( hASA_file, filename.str(), std::ios::out ) ) {
			hASA_file << chemical::name_from_aa( (chemical::AA)ii ) << "nbs21-24" << std::endl;
			for ( Size jj=1; jj <= hASA_by_aa_type_nbs_21_to_24[ ii ].size(); ++jj ) {
				hASA_file << hASA_by_aa_type_nbs_21_to_24[ ii ][ jj ] << std::endl;
			}
			hASA_file.close();
		} else {
			TR << "Could not open filename " << filename.str() << std::endl;
		}
		filename.str("");
	}

} // end calculate_hASA_by_type_and_nbcount


/// @details
/// Calculates patch areas for all surface positions. Surface positions are those with lte 16/20 nbs.
/// Uses all neighbor types, not just hydrophobic neighbors, or just surface-exposed residues.
/// Uses the average hASA residue values determined from lots of residue observations, not exact SASA values.
/// Does still include a check to make sure neighbors are surface-exposed. (If we didn't have this check, then the
/// patch areas would be significantly off.)
/// Splits up the patch areas by number of nbs that each surface position has for a finer grained score.
///
void
calculate_total_hASA_within_distance_avgresiduevalues_allnbs_conditionalonnumnbs( std::vector< FileName > & pdb_file_names ) {

	// create a score function using the standard packer weights
	// this will be used only to get a tenA neighbor graph
		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// keep track of the amount of hydrophobic ASA every surface-exposed residue has
	utility::vector1< utility::vector1< Real > > hASA_within_10A( 20 );

	// read in the per-residue hASA means
	// they'll be accessible via this map like: map[ ALA-lte10 ], map[ GLN-lte13 ], map[ TYR-lte16 ]
	std::map< std::string, Real > res_to_average_hASA;

	utility::io::izstream residue_hASA_ifstream;
	basic::database::open( residue_hASA_ifstream, "scoring/score_functions/SurfacePotential/average_hASA_by_res_and_neighbor.txt" );

	std::string restype;
	Real lte10_asa = 0.0;
	Real lte13_asa = 0.0;
	Real lte16_asa = 0.0;
	Real lte20_asa = 0.0;
	Real lte24_asa = 0.0;

	while ( !residue_hASA_ifstream.eof() ) {
		residue_hASA_ifstream >> restype >> lte10_asa >> lte13_asa >> lte16_asa >> lte20_asa >> lte24_asa;

		res_to_average_hASA[ restype + "-lte10" ] = lte10_asa;
		res_to_average_hASA[ restype + "-lte13" ] = lte13_asa;
		res_to_average_hASA[ restype + "-lte16" ] = lte16_asa;
		res_to_average_hASA[ restype + "-lte20" ] = lte20_asa;
		res_to_average_hASA[ restype + "-lte24" ] = lte24_asa;
	}

	// iterate through all the structures - do something to them
	for ( std::vector< FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb ) {

		TR << "Processing " << *pdb << "... " << std::endl;
		if ( !utility::file::file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
		if ( pose.size() <= 4 ) {
			std::cerr << "PDB not big enough. Quitting.";
			exit(1);
		}

		scorefxn->score( pose );

		Real patchArea = 0.0;
		Size countNeighbors = 0;

		// Now go through every residue (except the termini) and if it's on the surface, add that residues's hASA (based on its
		// neighbor count) and all neighboring residues hASA's to the total.  Save that final sum in the hASA_within_10A vector.
		// Do we want to include only the surface-exposed neighbors or should I include all neighbors since the average
		// hydrophobic ASA values include a neighbor dependence to them.

		for ( Size ii=4; ii <= pose.size()-3; ++ii ) {

			if ( ! pose.residue( ii ).is_protein() )
				continue;

			// our definition of surface residue is that the residue has fewer than 16 (or 20) neighbors
			// only calculate the hpatch data for residues with fewer than that # of nbs.
			// but neighboring residues which will contribute to this residues patch area should be allowed to have up to 24 nbs
			// before they contribute no hASA.
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
			countNeighbors = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
			if ( countNeighbors > SURFACE_EXPOSED_CUTOFF )
				continue;

			// passed the surface-exposed check...
			// reset the area size for every surface-exposed residue
			patchArea = 0.0;

			// add in the hASA contribution even if the residue is polar
			std::string key = get_map_key( pose.residue(ii).name3(), countNeighbors );
			if ( res_to_average_hASA.find( key ) != res_to_average_hASA.end() ) {
				patchArea += res_to_average_hASA.find( key )->second;
				#ifdef FILE_DEBUG
					TR << pose.residue( ii ).name3() << " " << ii << ": nbs: " << countNeighbors << ", hp_ASA: " << res_to_average_hASA.find( key )->second
						<< "; pdb: " << pose.pdb_info()->name() << std::endl;
				#endif
			} else {
				TR << "Residue " << pose.residue( ii ).name3() << " " << ii << ", nbs: " << countNeighbors << ", no hASA info; residue key: " << key << ", pdb: "
					<< pose.pdb_info()->name() << std::endl;
				exit(1);
			}

			// add in every neighbors hASA too
			Size nbresidue_num_nbs = 0; // create another variable to hold variable count because we want to use the one above for array index below
			#ifdef FILE_DEBUG
				TR << "Neighbors of residue " << pose.residue(ii).name3() << " " << ii << " ( " << res_to_average_hASA.find( key )->second << " ) include " << std::endl;
			#endif
			for ( utility::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node(ii)->const_edge_list_begin(),
				eli_end = tenA_neighbor_graph.get_node(ii)->const_edge_list_end(); eli != eli_end; ++eli ) {

				// get the other node for this edge, so pass in the res1 node to this method
				// save the value to simplify code ahead
				Size jj = (*eli)->get_other_ind(ii);

				if ( ! pose.residue(jj).is_protein() ) { continue; }

				#ifdef FILE_DEBUG
					TR << pose.residue(jj).name3() << " " << jj;
				#endif

				nbresidue_num_nbs = tenA_neighbor_graph.get_node(jj)->num_neighbors_counting_self();
				if ( nbresidue_num_nbs > BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
					#ifdef FILE_DEBUG
						TR << std::endl;
					#endif
					continue; // accomplishes the same as adding 0.00 for this residue/nb
				} else {
					key = get_map_key( pose.residue(jj).name3(), nbresidue_num_nbs );
					if ( res_to_average_hASA.find( key ) != res_to_average_hASA.end() ) {
						#ifdef FILE_DEBUG
							TR << ": nbs: " << nbresidue_num_nbs << ", hASA: " << res_to_average_hASA.find( key )->second << "; pdb: " << pose.pdb_info()->name();
						#endif
						patchArea += res_to_average_hASA.find( key )->second;
					} else {
						TR << pose.residue(jj).name3() << " " << jj << ", nbs: " << nbresidue_num_nbs << ", no hASA info; residue key: " << key << std::endl;
						exit(1);
					}
				}
				#ifdef FILE_DEBUG
					TR << std::endl;
				#endif
			}
			#ifdef FILE_DEBUG
				TR << "\n" << std::endl;
			#endif

			// now that we know how much exposed hASA this residue has, save it somewhere
			hASA_within_10A[ countNeighbors ].push_back( patchArea );

		} // end loop over all residues

	} // end for loop over multiple input pdb files

	#ifdef FILE_DEBUG
		TR << "Done looping over PDB files. Outputting patch files." << std::endl;
	#endif

	// output the data into separate files for each residue type
	std::stringstream filename;
	std::ofstream patches_file;

	for ( Size ii=1; ii <= hASA_within_10A.size(); ++ii ) {
		filename << "hASA_patch_sizes.avgresiduevalues.allnbs.txt.nb" << I(2,ii);
		if ( trytry_ofstream_open( patches_file, filename.str(), std::ios::out ) ) {
			patches_file << "nb" << ii << std::endl;
			for ( Size jj=1; jj <= hASA_within_10A[ ii ].size(); ++jj ) {
				patches_file << hASA_within_10A[ ii ][ jj ] << std::endl;
			}
			patches_file.close();
			filename.str("");
		} else {
			TR << "Could not open filename " << filename << std::endl;
		}
	}

} // end calculate_total_hASA_within_distance_avgresiduevalues_allnbs_conditionalonnumnbs


/// @brief main function
int
main( int argc, char* argv[] )
{

	try {

	clock_t starttime = clock();
	//basic::prof_reset();

	// options, random initialization
	devel::init( argc, argv );

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//ronj concatenate -s and -l flags together to get total list of PDB files
	//ronj  The advantage of parsing -s and -l separately is that users can specify a list and a single structure on the
	//ronj command line.
	std::vector< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)

	std::vector< FileName > list_file_names;
	if ( option[ in::file::l ].active() ) {
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)

		//ronj for each file input with the -l switch...
		for ( std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i ) {
			std::string listfilename( i->name() );
			std::ifstream data( listfilename.c_str() );
			//ronj try to open that particular file first...
			if ( !data.good() ) {
				utility_exit_with_message( "Unable to open file: " + listfilename + '\n' );
			}
			std::string line;
			//ronj then read all the lines in that file until there are no more
			while( getline(data, line) ) {
				pdb_file_names.push_back( FileName(line) );
			}
			data.close();
		}
	}

	if ( pdb_file_names.size() == 0 ) {
		utility_exit_with_message_status("No files given: Use either -file:s or -file:l to designate a single pdb or a list of pdbs", 1);
	}

	TR << "Constructed list of pdb files." << std::endl;

	//find_hppatches_nb_graph( pdb_file_names );
	//find_hppatches_distance( pdb_file_names );
	//calculate_percent_hydrophobic_stats( pdb_file_names );
	//calculate_percent_hydrophobic_distribution( pdb_file_names );
	//calculate_hydrophobic_accessible_surface_area( pdb_file_names );
	//calculate_total_hASA_within_distance_exact_hASA_values( pdb_file_names );
	//calculate_total_hASA_within_distance_exact_hASA_values_allnbs( pdb_file_names );

	//calculate_total_hASA_within_distance_miniSASAvalues_allnbs_exposedornot_conditionalonnumnbs( pdb_file_names );
	//calculate_total_hASA_within_distance_exact_hASA_values_allnbs_exposedornot_conditionalonnumnbs( pdb_file_names );

	//calculate_total_hASA_within_distance_avgresiduevalues_allnbs_conditionalonnumnbs( pdb_file_names );
	//calculate_hASA_by_type_and_exposure( pdb_file_names );
	//calculate_hASA_by_type_and_attractiveE( pdb_file_names );

	//calculate_hydrophobic_accessible_surface_area( pdb_file_names );

	//calculate_hASA_by_type_and_nbcount( pdb_file_names );
	calculate_total_hASA_within_distance_avgresiduevalues_allnbs_conditionalonnumnbs( pdb_file_names );

	//basic::prof_show();
	clock_t stoptime = clock();
	TR << "Whole run took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

