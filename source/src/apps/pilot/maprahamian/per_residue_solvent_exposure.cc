// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// /// @file per_residue_solvent_exposure.cc
// /// @brief Determines the per residue solvent exposure using one of the following methods:
//       1) CENTROID
//        a) sphere method ( -centroid_version -sphere_method)
//        b) cone method ( -centroid_version -cone_method)
//       2) FULL ATOM
//        a) closest atom to target CA counted as neighbor
//         i) sphere method ( -neighbor_closest_atom -sphere_method)
//         ii) cone method ( -neighbor_closest_atom -cone_method)
//        b) CB atom counted as neighbor
//         i) sphere method ( -sphere_method)
//         ii) cone method ( -cone_method)
//        c) number of neighboring atoms
//         i) sphere method ( -atom_neighbor_count -sphere_method)
//         ii) cone method ( -atom_neighbor_count -cone_method)
//       Additional flags:
//        -dist_midpoint (distance logistic function midpoint)
//        -dist_steepness (exponential falloff of distance function)
//        -angle_midpoint (angle logistic function midpoint)
//        -angle_steepness (exponential falloff of angle function)
//        -cone_angle (angle cutoff for the cone)
// /// @author Melanie Aprahamian (aprahamian.4@osu.edu)

#include <iostream>
#include <string>
#include <protocols/jd2/util.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <numeric/NumericTraits.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;

OPT_KEY( Real, dist_midpoint )
OPT_KEY( Real, dist_steepness )
OPT_KEY( Real, cone_angle )
OPT_KEY( Real, angle_midpoint )
OPT_KEY( Real, angle_steepness )
OPT_KEY( Boolean, sphere_method )
OPT_KEY( Boolean, cone_method )
OPT_KEY( Boolean, centroid_version )
OPT_KEY( Boolean, neighbor_closest_atom )
OPT_KEY( Boolean, atom_neighbor_count )

static basic::Tracer TR( "apps.pilot.maprahamian.per_residue_solvent_exposure" );

// Main program
int main( int argc, char * argv [] ){
	try {

		NEW_OPT( dist_midpoint, "midpoint of distance falloff", 9.0);
		NEW_OPT( dist_steepness, "exponent factor for distance falloff", 1.0);
		NEW_OPT( cone_angle, "angle cutoff for cone", numeric::NumericTraits< float >::pi_2());
		NEW_OPT( angle_midpoint, "midpoint of angle falloff", numeric::NumericTraits< float >::pi_2()/2.0);
		NEW_OPT( angle_steepness, "exponent factor for angle falloff", numeric::NumericTraits< float >::pi()*2.0);
		NEW_OPT( sphere_method, "calculate neighbor count based upon SPHERE method", false);
		NEW_OPT( cone_method, "calculate neighbor count based upon the CONE method", false);
		NEW_OPT( centroid_version, "centroid versions of the neighbor counts", false);
		NEW_OPT( neighbor_closest_atom, "for FA neighbor count, use closest atom per residue as neighbor", false);
		NEW_OPT( atom_neighbor_count, "atom neighbor count", false);

		// Initialize Rosetta
		devel::init( argc, argv );

		// Namespaces
		using namespace core;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace pose;

		// check if input PDB file is provided from command line with option -s
		if ( !option[in::file::s].user() ) {
			// exit if no PDB file is found
			utility_exit_with_message("Input PDB file not found, please use the -s flag");
		}

		core::Real distance_fn_midpoint = option[ dist_midpoint ];
		core::Real distance_fn_steepness = option[ dist_steepness ];
		core::Real cone_angle_cutoff = option[ cone_angle ];
		core::Real angle_fn_midpoint = option[ angle_midpoint ];
		core::Real angle_fn_steepness = option[ angle_steepness ];
		bool sphere = option[ sphere_method ];
		bool cone = option[ cone_method ];
		bool centroid_pose = option[ centroid_version ];
		bool closest_atom = option[ neighbor_closest_atom ];
		bool atom_neighbors = option[ atom_neighbor_count ];

		pose::Pose pose;
		std::string pdb_file = option[in::file::s]()[1];

		//////////////////////
		// CENTROID versions
		if ( centroid_pose == true ) {
			core::import_pose::centroid_pose_from_pdb( pose, pdb_file, core::import_pose::PDB_file);

			// SPHERE method
			if ( sphere == true ) {
				core::Size const num_residues (pose.total_residue());
				for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
					core::Real neighbor_count (0.0);
					for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
						if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
							core::Real const distance = pose.residue(res_count_neighbor).xyz("CEN").distance(pose.residue(res_count_target).xyz("CEN"));
							neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
						}
					}
					TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
				}
			} // SPHERE

			// CONE method
			if ( cone == true ) {
				core::Size const num_residues (pose.total_residue());
				for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
					core::Real neighbor_count (0.0);
					numeric::xyzVector<core::Real> CA_CEN_vector = pose.residue(res_count_target).xyz("CEN") - pose.residue(res_count_target).xyz("CA");
					core::Real const distance_internal = pose.residue(res_count_target).xyz("CEN").distance(pose.residue(res_count_target).xyz("CA"));
					for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
						if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
							numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz("CEN") - pose.residue(res_count_target).xyz("CA");
							core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz("CEN").distance(pose.residue(res_count_target).xyz("CA"));
							numeric::xyzVector<core::Real> norm_CA_CEN_vector = CA_CEN_vector/distance_internal;
							numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
							core::Real angle = std::acos(norm_CA_CEN_vector.dot(norm_neighbor_vector));
							if ( angle > cone_angle_cutoff ) {
								angle = 0.0;
							}
							neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
						}
					}
					TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
				}
			} // CONE
		} // CENTROID

		/////////////////////////////
		//
		// FULL ATOM versions
		if ( centroid_pose == false ) {
			core::import_pose::pose_from_file(pose,pdb_file);

			// RESIDUE NEIGHBORS version
			if ( atom_neighbors == false ) {
				// SPHERE method
				if ( sphere == true ) {
					core::Size const num_residues (pose.total_residue());
					for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
						core::Real neighbor_count (0.0);
						std::string CB_atom ("CB");
						if ( pose.residue(res_count_target).type().name1() == 'G' ) {
							CB_atom = "1HA";
						}
						for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
							if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
								std::string neighbor_atom ("CB");
								if ( pose.residue(res_count_neighbor).type().name1() == 'G' ) {
									neighbor_atom = "1HA";
								}
								if ( closest_atom == true ) {
									std::string neighbor_atom (pose.residue(res_count_neighbor).atom_name(1));
									for ( core::Size atom_count = 2; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
										if ( pose.residue(res_count_neighbor).xyz(atom_count).distance(pose.residue(res_count_target).xyz(CB_atom)) < pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(CB_atom)) ) {
											neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
										}
									}
								}
								core::Real const distance = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(CB_atom));
								neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
							}
						}
						TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
					}
				} // SPHERE

				// CONE method
				if ( cone == true ) {
					core::Size const num_residues (pose.total_residue());
					for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
						core::Real neighbor_count (0.0);
						std::string CB_atom ("CB");
						if ( pose.residue(res_count_target).type().name1() == 'G' ) {
							CB_atom = "1HA";
						}
						numeric::xyzVector<core::Real> CA_CB_vector = pose.residue(res_count_target).xyz(CB_atom) - pose.residue(res_count_target).xyz("CA");
						core::Real const distance_internal = pose.residue(res_count_target).xyz(CB_atom).distance(pose.residue(res_count_target).xyz("CA"));
						for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
							if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
								std::string neighbor_atom ("CB");
								if ( pose.residue(res_count_neighbor).type().name1() == 'G' ) {
									neighbor_atom = "1HA";
								}
								if ( closest_atom == true ) {
									std::string neighbor_atom (pose.residue(res_count_neighbor).atom_name(1));
									for ( core::Size atom_count = 2; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
										if ( pose.residue(res_count_neighbor).xyz(atom_count).distance(pose.residue(res_count_target).xyz(CB_atom)) < pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(CB_atom)) ) {
											neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
										}
									}
								}
								numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz(neighbor_atom) - pose.residue(res_count_target).xyz("CA");
								core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz("CA"));
								numeric::xyzVector<core::Real> norm_CA_CB_vector = CA_CB_vector/distance_internal;
								numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
								core::Real angle = std::acos(norm_CA_CB_vector.dot(norm_neighbor_vector));
								if ( angle > cone_angle_cutoff ) {
									angle = 0.0;
								}
								neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
							}
						}
						TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
					}
				} // CONE
			} // RESIDUE NEIGHBORS

			// ATOM NEIGHBORS version
			if ( atom_neighbors == true ) {
				// SPHERE method
				if ( sphere == true ) {
					core::Size const num_residues (pose.total_residue());
					for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
						core::Real neighbor_count (0.0);
						std::string CB_atom ("CB");
						if ( pose.residue(res_count_target).type().name1() == 'G' ) {
							CB_atom = "1HA";
						}
						for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
							if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
								for ( core::Size atom_count = 1; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
									std::string neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
									core::Real const distance = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(CB_atom));
									neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
								}
							}
						}
						TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
					}
				} // SPHERE

				// CONE method
				if ( cone == true ) {
					core::Size const num_residues (pose.total_residue());
					for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
						core::Real neighbor_count (0.0);
						std::string CB_atom ("CB");
						if ( pose.residue(res_count_target).type().name1() == 'G' ) {
							CB_atom = "1HA";
						}
						numeric::xyzVector<core::Real> CA_CB_vector = pose.residue(res_count_target).xyz(CB_atom) - pose.residue(res_count_target).xyz("CA");
						core::Real const distance_internal = pose.residue(res_count_target).xyz(CB_atom).distance(pose.residue(res_count_target).xyz("CA"));
						for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
							if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
								for ( core::Size atom_count = 1; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
									std::string neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
									numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz(neighbor_atom) - pose.residue(res_count_target).xyz("CA");
									core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz("CA"));
									numeric::xyzVector<core::Real> norm_CA_CB_vector = CA_CB_vector/distance_internal;
									numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
									core::Real angle = std::acos(norm_CA_CB_vector.dot(norm_neighbor_vector));
									if ( angle > cone_angle_cutoff ) {
										angle = 0.0;
									}
									neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
								}
							}
						}
						TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
					}
				} // CONE
			} // ATOM NEIGHBORS
		} // FULL ATOM
	}
catch ( utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
} //main


