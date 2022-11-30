// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief aginsparg, ipatel, sthyme


#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/motifs/Motif.hh>
#include <core/chemical/residue_io.hh>
#include <protocols/motifs/MotifHit.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/motifs/BuildPosition.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <cmath>
#include <core/pose/PDBInfo.hh>

#include <core/chemical/GlobalResidueTypeSet.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
using namespace protocols::dna;

#include <basic/prof.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR( "apps.pilot.motif_dna_packer_design" );

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
using utility::string_split;

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

//this may solve atom clash issue
#include <core/pose/xyzStripeHashPose.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>

#include <numeric/xyzVector.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/ligand_scores.hh>

#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <ObjexxFCL/FArray1D.hh>

// Time profiling header
#include <time.h>


//for data type debugging
#include <typeinfo>

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;

////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		static basic::Tracer ms_tr( "Code tracing", basic::t_info );

		//set up clock for recording program runtime
		clock_t program_start = clock();

		using namespace options;
		using namespace OptionKeys;

		devel::init( argc, argv ); // reading options--name should be more descriptive

		basic::prof_reset();

		ms_tr << "Making Pose of input pdb file(s)" << std::endl;

		typedef vector1< utility::file::FileName > Filenames;
		Filenames pdbnames;

		//borrowing code from motif_ligand_packer_design.cc to take in a single or multiple pdb files (using -l or -s flags), and be able to iterate discovery upon all files
		if ( option[ in::file::l ].user() ) {
			Filenames listnames( option[ in::file::l ]().vector() );
			for ( Filenames::const_iterator filename( listnames.begin() );
					filename != listnames.end(); ++filename ) {
				std::ifstream list( (*filename).name().c_str() );
				while ( list ) {
					std::string pdbname;
					list >> pdbname;
					pdbnames.push_back( pdbname );
				}
			}

		} else if ( option[ in::file::s ].user() ) {
			pdbnames = option[ in::file::s ]().vector();

		} else {
			std::cerr << "No files given: Use either -file:s or -file:l "
				<< "to designate a single pdb or a list of pdbs"
				<< std::endl;
		}

		//run space_filling for each pdb
		for ( Filenames::const_iterator filename( pdbnames.begin() );
				filename != pdbnames.end(); ++filename ) {
			if ( !utility::file::file_exists( *filename ) ) {
				continue;
			}

			pose::Pose pose_pointer_helper;

			//ms_tr << "Before import pose" << std::endl;
			core::import_pose::pose_from_file( pose_pointer_helper, *filename , core::import_pose::PDB_file);
			pose::PoseOP pose(new pose::Pose(pose_pointer_helper));

			//ms_tr << pose_pointer_helper.size() << ", " << pose->size() << std::endl;
			//ms_tr << "After import pose" << std::endl;
			std::string pdbprefix( string_split( string_split( *filename, '/' ).back(), '.' ).front() );



			//currently deeclaring target_conformers_map here; will need to define earlier in a class declaration when converting from app to protocol

			//run through all atoms to derive a range of dimensions to contain the protein in a 3D  space
			//since we can't have negative indices, we need to normalize the coordinate values so that everything is positive
			//derive constant values based  on the most negative values in each dimension, and then add that constant to all coordinates

			int smallest_x = 1;
			int smallest_y = 1;
			int smallest_z = 1;

			int largest_x = 1;
			int largest_y = 1;
			int largest_z = 1;

			//variable to  hold the linear scale by which we will increase the resolution of identification
			//default is 1. Probably want to set to 3-6
			int resolution_increase_factor = option[ OptionKeys::motifs::resolution_scale_factor ];

			//create a list of coordinates of each atom to hold and work with to fill the protein_representation_matrix
			//can't seem to make a vector of xyzVector objects, so will need to just make a custome 2D vector  to  hold the data
			//unlike clash detection, going to use floats as well for increaserd precision
			utility::vector1<numeric::xyzVector<int>> atom_coordinates;


			utility::vector1<utility::vector1<core::Real>> atom_coordinates_float_and_lj_radius;

			//determine largest and smallest x,y,z  values to determine dimensions of matrix
			for ( core::Size res_num = 1; res_num <= pose->size(); ++res_num ) {
				for ( core::Size atom_num = 1; atom_num <= pose->residue(res_num).natoms(); ++atom_num ) {
					//get the x,y,z data of the atom
					numeric::xyzVector<int> atom_xyz = pose->residue(res_num).xyz(atom_num);

					utility::vector1<core::Real> atom_xyz_float_with_lj_radius;
					atom_xyz_float_with_lj_radius.push_back(atom_xyz.x());
					atom_xyz_float_with_lj_radius.push_back(atom_xyz.y());
					atom_xyz_float_with_lj_radius.push_back(atom_xyz.z());
					atom_xyz_float_with_lj_radius.push_back(pose->residue(res_num).atom_type(atom_num).lj_radius());

					atom_coordinates_float_and_lj_radius.push_back(atom_xyz_float_with_lj_radius);


					//convert the coordinates to integers
					atom_xyz.x() = static_cast<int>(atom_xyz.x());
					atom_xyz.y() = static_cast<int>(atom_xyz.y());
					atom_xyz.z() = static_cast<int>(atom_xyz.z());

					//determine if any of the values  are the smallest
					if ( smallest_x > atom_xyz.x() ) {
						smallest_x = atom_xyz.x();
					}
					if ( smallest_y > atom_xyz.y() ) {
						smallest_y = atom_xyz.y();
					}
					if ( smallest_z > atom_xyz.z() ) {
						smallest_z = atom_xyz.z();
					}

					//determine if any  are the largest
					if ( largest_x < atom_xyz.x() ) {
						largest_x = atom_xyz.x();
					}
					if ( largest_y < atom_xyz.y() ) {
						largest_y = atom_xyz.y();
					}
					if ( largest_z < atom_xyz.z() ) {
						largest_z = atom_xyz.z();
					}

					atom_coordinates.push_back(atom_xyz);

				}
			}

			//take negative values of the smallest values and then add 1 to derive the constants
			core::Size x_shift = (smallest_x * -1) + 1;
			core::Size y_shift = (smallest_y * -1) + 1;
			core::Size z_shift = (smallest_z * -1) + 1;

			//apply shift values to largest to get boundaries
			core::Size x_bound  = x_shift + largest_x;
			core::Size y_bound  = y_shift + largest_y;
			core::Size z_bound  = z_shift + largest_z;

			//int x_bound_int = x_bound * resolution_increase_factor;
			//int y_bound_int = y_bound * resolution_increase_factor;
			//int z_bound_int = z_bound * resolution_increase_factor;

			//apply constant shifts to all coordinates
			//apply resolution factor to all atoms
			for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates.size(); ++xyzVec ) {
				atom_coordinates[xyzVec].x() += x_shift;
				atom_coordinates[xyzVec].y() += y_shift;
				atom_coordinates[xyzVec].z() += z_shift;

				atom_coordinates[xyzVec].x() *= resolution_increase_factor;
				atom_coordinates[xyzVec].y() *= resolution_increase_factor;
				atom_coordinates[xyzVec].z() *= resolution_increase_factor;

				atom_coordinates_float_and_lj_radius[xyzVec][1] += x_shift;
				atom_coordinates_float_and_lj_radius[xyzVec][2] += y_shift;
				atom_coordinates_float_and_lj_radius[xyzVec][3] += z_shift;

				atom_coordinates_float_and_lj_radius[xyzVec][1] *= resolution_increase_factor;
				atom_coordinates_float_and_lj_radius[xyzVec][2] *= resolution_increase_factor;
				atom_coordinates_float_and_lj_radius[xyzVec][3] *= resolution_increase_factor;
				atom_coordinates_float_and_lj_radius[xyzVec][4] *= resolution_increase_factor;
			}

			//create 3D matrix to roughly represent 3D coordinate space of protein
			//bad syntax, may just have to do  an iterative  fill
			//utility::vector1<utility::vector1<utility::vector1<bool>>> protein_representation_matrix(x_bound, (y_bound, (z_bound, false)));

			x_bound *= resolution_increase_factor;
			y_bound *= resolution_increase_factor;
			z_bound *= resolution_increase_factor;

			ms_tr << "Creating protein clash coordinate matrix. Dimensions of matrix are " << x_bound << "," << y_bound << "," << z_bound << std::endl;

			utility::vector1<utility::vector1<utility::vector1<bool>>> protein_representation_matrix;

			for ( core::Size x = 1; x <= x_bound; ++x ) {
				//make a 2D  matrix
				utility::vector1<utility::vector1<bool>> sub_matrix;

				for ( core::Size y = 1; y <= y_bound; ++y ) {

					//make a 1D matrix, seed with false values
					utility::vector1<bool> sub_sub_matrix(z_bound,  false);
					//push 1D  matrix into 2D
					sub_matrix.push_back(sub_sub_matrix);

				}
				//push a 2D  matrix into the 3D matrix
				protein_representation_matrix.push_back(sub_matrix);
			}

			//seed the matrix with approximate coordinates of each atom
			//use the vector that has the LJ radii to fill out cells
			for ( core::Size xyzVec = 1; xyzVec <= atom_coordinates_float_and_lj_radius.size(); ++xyzVec ) {
				//iterate through all cells in the matrix and determine if the sphere projected by the atom center and LJ radius hits this cell
				//there may be a more efficient way to do this...
				//I think there is! instead of iterating the whole volume of the matrix, only iterate about a range bound by a cube with with side length = 2 x LJ radius
				//adjust min and max to ensure we don't run off of the matrix

				core::Size x_min = atom_coordinates_float_and_lj_radius[xyzVec][1] - atom_coordinates_float_and_lj_radius[xyzVec][4];
				core::Size x_max = atom_coordinates_float_and_lj_radius[xyzVec][1] + atom_coordinates_float_and_lj_radius[xyzVec][4];
				if ( x_min < 1 ) {
					x_min = 1;
				}
				if ( x_max > x_bound ) {
					x_max = x_bound;
				}
				core::Size y_min = atom_coordinates_float_and_lj_radius[xyzVec][2] - atom_coordinates_float_and_lj_radius[xyzVec][4];
				core::Size y_max = atom_coordinates_float_and_lj_radius[xyzVec][2] + atom_coordinates_float_and_lj_radius[xyzVec][4];
				if ( y_min < 1 ) {
					y_min = 1;
				}
				if ( y_max > y_bound ) {
					y_max = y_bound;
				}
				core::Size z_min = atom_coordinates_float_and_lj_radius[xyzVec][3] - atom_coordinates_float_and_lj_radius[xyzVec][4];
				core::Size z_max = atom_coordinates_float_and_lj_radius[xyzVec][3] + atom_coordinates_float_and_lj_radius[xyzVec][4];
				if ( z_min < 1 ) {
					z_min = 1;
				}
				if ( z_max > y_bound ) {
					z_max = z_bound;
				}

				for ( core::Size x = x_min; x <= x_max; ++x ) {
					for ( core::Size y = y_min; y <= y_max; ++y ) {
						for ( core::Size z = z_min; z <= z_max; ++z ) {
							//std::cout << atom_coordinates_float_and_lj_radius.size() << "," << xyzVec << "," << x << "," << y << "," << z << std::endl;
							//std::cout << x_min << "," << x_max << "," << x_bound << std::endl;
							//std::cout << y_min << "," << y_max << "," << y_bound << std::endl;
							//std::cout << z_min << "," << z_max << "," << z_bound << std::endl;

							//use distance formula to figure out if cell x,y,z is within the sphere projected by the atom point about its LJ radius
							//get distance between x,y,z and the atom point
							//std::cout << "Calculating distance" << std::endl;
							core::Real atom_cell_distance = sqrt((x - atom_coordinates_float_and_lj_radius[xyzVec][1]) + (y - atom_coordinates_float_and_lj_radius[xyzVec][2]) + (z - atom_coordinates_float_and_lj_radius[xyzVec][3]));
							//std::cout << "Distance calculated" << std::endl;
							//if distance is less than the radius, then the point is occupied
							if ( atom_cell_distance < (x - atom_coordinates_float_and_lj_radius[xyzVec][4]) ) {
								//std::cout << "Occupied" << std::endl;
								protein_representation_matrix[x][y][z] = true;
								//std::cout << "Cell marked" << std::endl;
							}
						}
					}
				}
			}

			//system should be processed. Get a count of occupied vs unoccupied cells
			core::Real total_cells = x_bound * y_bound * z_bound;

			core::Real occupied_cell_count = 0;
			core::Real unoccupied_cell_count = 0;



			for ( core::Size x = 1; x <= x_bound; ++x ) {
				for ( core::Size y = 1; y <= y_bound; ++y ) {
					for ( core::Size z = 1; z <= z_bound; ++z ) {
						if ( protein_representation_matrix[x][y][z] == true ) {
							++occupied_cell_count;
						} else {
							++unoccupied_cell_count;
						}
					}
				}
			}

			core::Real occupied_ratio = occupied_cell_count / total_cells;

			std::cout << "Total: " << total_cells << std::endl;
			std::cout << "Occupied: " << occupied_cell_count << std::endl;
			std::cout << "Unoccupied: " << unoccupied_cell_count << std::endl;
			std::cout << "Occupied-Total Ratio: " << occupied_ratio << std::endl;

			//investigate a specified sub_region of the system (only if values are given in the arguments)
			if ( option[ OptionKeys::motifs::x_min_holes ].user() && option[ OptionKeys::motifs::x_max_holes ].user() && option[ OptionKeys::motifs::y_min_holes ].user() && option[ OptionKeys::motifs::y_max_holes ].user() && option[ OptionKeys::motifs::z_min_holes ].user() && option[ OptionKeys::motifs::z_max_holes ].user() ) {
				std::cout << "Options entered to get info on a specific boundary within the PDB. Using initially-inputted boundaries: " << std::endl;

				//get the xyz minima and maxima of the region in the un-motified format
				core::Real x_min_holes = option[ OptionKeys::motifs::x_min_holes ];
				core::Real x_max_holes = option[ OptionKeys::motifs::x_max_holes ];
				core::Real y_min_holes = option[ OptionKeys::motifs::y_min_holes ];
				core::Real y_max_holes = option[ OptionKeys::motifs::y_max_holes ];
				core::Real z_min_holes = option[ OptionKeys::motifs::z_min_holes ];
				core::Real z_max_holes = option[ OptionKeys::motifs::z_max_holes ];

				std::cout << x_min_holes << "," << x_max_holes << std::endl;
				std::cout << y_min_holes << "," << y_max_holes << std::endl;
				std::cout << z_min_holes << "," << z_max_holes << std::endl;

				//apply shift to these boundaries and then the resolution factor
				x_min_holes += x_shift;
				x_max_holes += x_shift;
				y_min_holes += y_shift;
				y_max_holes += y_shift;
				z_min_holes += z_shift;
				z_max_holes += z_shift;

				x_min_holes *= resolution_increase_factor;
				x_max_holes *= resolution_increase_factor;
				y_min_holes *= resolution_increase_factor;
				y_max_holes *= resolution_increase_factor;
				z_min_holes *= resolution_increase_factor;
				z_max_holes *= resolution_increase_factor;

				std::cout << "After shifts and resolution factor:" << std::endl;
				std::cout << x_min_holes << "," << x_max_holes << std::endl;
				std::cout << y_min_holes << "," << y_max_holes << std::endl;
				std::cout << z_min_holes << "," << z_max_holes << std::endl;

				//cast values to integers for matrix accession
				//ensure that values are within 1 and respective boundary
				core::Size x_min_holes_int = x_min_holes;
				if ( x_min_holes_int < 1 ) {
					x_min_holes = 1;
				}
				if ( x_min_holes_int > x_bound ) {
					x_min_holes = x_bound;
				}
				core::Size x_max_holes_int = x_max_holes;
				if ( x_max_holes_int < 1 ) {
					x_max_holes = 1;
				}
				if ( x_max_holes_int > x_bound ) {
					x_max_holes = x_bound;
				}
				core::Size y_min_holes_int = y_min_holes;
				if ( y_min_holes_int < 1 ) {
					y_min_holes = 1;
				}
				if ( y_min_holes_int > y_bound ) {
					y_min_holes = y_bound;
				}
				core::Size y_max_holes_int = y_max_holes;
				if ( y_max_holes_int < 1 ) {
					y_max_holes = 1;
				}
				if ( y_max_holes_int > y_bound ) {
					y_max_holes = y_bound;
				}
				core::Size z_min_holes_int = z_min_holes;
				if ( z_min_holes_int < 1 ) {
					z_min_holes = 1;
				}
				if ( z_min_holes_int > z_bound ) {
					z_min_holes = z_bound;
				}
				core::Size z_max_holes_int = z_max_holes;
				if ( z_max_holes_int < 1 ) {
					z_max_holes = 1;
				}
				if ( z_max_holes_int > z_bound ) {
					z_max_holes = z_bound;
				}

				std::cout << "Integers used and accounting to ensure we don't go out of bounds" << std::endl;
				std::cout << x_min_holes_int << "," << x_max_holes_int << std::endl;
				std::cout << y_min_holes_int << "," << y_max_holes_int << std::endl;
				std::cout << z_min_holes_int << "," << z_max_holes_int << std::endl;

				//look through portion of the matrix and get volume of occupied vs unoccupied

				core::Real subsection_occupied = 0;
				core::Real subsection_unoccupied = 0;
				core::Real subsection_total = 0;

				for ( core::Size x = x_min_holes; x <= x_max_holes; ++x ) {
					for ( core::Size y = y_min_holes; y <= y_max_holes; ++y ) {
						for ( core::Size z = z_min_holes; z <= z_max_holes; ++z ) {
							++subsection_total;
							if ( protein_representation_matrix[x][y][z] == true ) {
								++subsection_occupied;
							} else {
								++subsection_unoccupied;
							}
						}
					}
				}

				std::cout << "Total: " << subsection_total << std::endl;
				std::cout << "Occupied: " << subsection_occupied << std::endl;
				std::cout << "Unoccupied: " << subsection_unoccupied << std::endl;

				core::Real subsection_ratio = subsection_occupied / subsection_total;

				std::cout << "Occupied-Total Ratio: " << subsection_ratio << std::endl;
			}



		}//end of iteration through all pdbs

		ms_tr << "SUCCESSFUL COMPLETION" << std::endl;

		clock_t program_total_time = clock() - program_start;

		std::cout << "Program completed in: " << ((float)program_total_time)/CLOCKS_PER_SEC << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
