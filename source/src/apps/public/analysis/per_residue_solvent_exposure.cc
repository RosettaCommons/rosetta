// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file per_residue_solvent_exposure.cc
/// @brief Determines the per residue solvent exposure using one of the following methods:
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
//       Additional flags:
//        -dist_midpoint (distance logistic function midpoint)
//        -dist_steepness (exponential falloff of distance function)
//        -angle_midpoint (angle logistic function midpoint)
//        -angle_steepness (exponential falloff of angle function)
//        -cone_angle (angle cutoff for the cone)
//        -bbamide (using the N-H vector rather the CA-CB vector)
/// @author Melanie Aprahamian (aprahamian.4@osu.edu) and Daniel Marzolf (marzolf.4@osu.edu)

#include <iostream>
#include <string>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/solvent_exposure.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


#include <numeric/NumericTraits.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;

OPT_KEY( Real, dist_midpoint )
OPT_KEY( Real, dist_steepness )
OPT_KEY( Real, angle_midpoint )
OPT_KEY( Real, angle_steepness )
OPT_KEY( StringVector, neighbor_method )
OPT_KEY( Boolean, centroid_version )
OPT_KEY( Boolean, neighbor_closest_atom )
OPT_KEY( Boolean, bbamide )


static basic::Tracer TR( "apps.public.analysis.per_residue_solvent_exposure" );


//Output Results; overwrites previous/existing file
void output_results( std::string const & outf ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile = "NeighborCount_default.out";
	if ( option[ out::file::o ].user() ) {
		outfile = option[ out::file::o ]();
	}
	utility::io::ozstream outy( outfile );
	outy << outf;
	TR << "Neighbor Count Finished, find results at " << outfile << std::endl;
}

//Neighbor Count Function for BB_Amide Cone
core::Real count_neighbors( std::string neighbor_oxygen, pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness, core::Real angle_fn_midpoint, core::Real angle_fn_steepness, core::Size res_count_target, std::string target_atom_vector_start,  core::Real const distance_internal, numeric::xyzVector<core::Real> target_vector, core::Size res_count_neighbor ) {
	numeric::xyzVector<core::Real> neighbor_O_vector = pose.residue(res_count_neighbor).xyz(neighbor_oxygen) - pose.residue(res_count_target).xyz(target_atom_vector_start);
	core::Real const distance_to_O = pose.residue(res_count_neighbor).xyz(neighbor_oxygen).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
	numeric::xyzVector<core::Real> norm_target_vector = target_vector/distance_internal;
	numeric::xyzVector<core::Real> norm_O_vector = neighbor_O_vector/distance_to_O;
	core::Real angleO = std::acos(norm_target_vector.dot(norm_O_vector));
	core::Real resi_nc = 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_O-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angleO-angle_fn_midpoint)));
	return resi_nc;
}


//Sphere FA Method
void sphere_FA_nc_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		std::string target_atom ("CA");
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				std::string neighbor_atom ("CB");
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_gly ) {
					neighbor_atom = "1HA";
				}
				core::Real const distance = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}

// Sphere FA Closest Atom Method
void sphere_FA_nc_anyatom_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		std::string target_atom ("CA");
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				std::string neighbor_atom (pose.residue(res_count_neighbor).atom_name(1));
				for ( core::Size atom_count = 2; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
					if ( pose.residue(res_count_neighbor).xyz(atom_count).distance(pose.residue(res_count_target).xyz(target_atom)) < pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom)) ) {
						neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
					}
				}
				core::Real const distance = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}

// Sphere Centroid Method
void sphere_cen_nc_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	std::string target ("CEN");
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				core::Real const distance = pose.residue(res_count_neighbor).xyz("CEN").distance(pose.residue(res_count_target).xyz(target));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}

// Cone FA Method CA to CB
void cone_FA_nc_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness, core::Real angle_fn_midpoint, core::Real angle_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		std::string target_atom_vector_start ("CA");
		std::string target_atom ("CB");
		if ( pose.residue(res_count_target).aa() == core::chemical::AA::aa_gly ) {
			target_atom = "1HA";
		}
		numeric::xyzVector<core::Real> target_vector = pose.residue(res_count_target).xyz(target_atom) - pose.residue(res_count_target).xyz(target_atom_vector_start);
		core::Real const distance_internal = pose.residue(res_count_target).xyz(target_atom).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				std::string neighbor_atom ("CB");
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_gly ) {
					neighbor_atom = "1HA";
				}
				numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz(neighbor_atom) - pose.residue(res_count_target).xyz(target_atom_vector_start);
				core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
				numeric::xyzVector<core::Real> norm_target_vector = target_vector/distance_internal;
				numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
				core::Real angle = std::acos(norm_target_vector.dot(norm_neighbor_vector));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}

void cone_FA_nc_calculator_amide ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness, core::Real angle_fn_midpoint, core::Real angle_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		std::string target_atom_vector_start ("N");
		std::string target_atom ("H");
		if ( res_count_target == 1 ) {
			target_atom = "2H";
		}
		if ( pose.residue(res_count_target).aa() == core::chemical::AA::aa_pro ) {
			target_atom = "CB";
		}
		numeric::xyzVector<core::Real> target_vector = pose.residue(res_count_target).xyz(target_atom) - pose.residue(res_count_target).xyz(target_atom_vector_start);
		core::Real const distance_internal = pose.residue(res_count_target).xyz(target_atom).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				std::string neighbor_O ("O");
				neighbor_count += count_neighbors(neighbor_O, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_asn ) {
					std::string neighbor_OD1 ("OD1") ;
					neighbor_count += count_neighbors(neighbor_OD1, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_tyr ) {
					std::string neighbor_OH ("OH") ;
					neighbor_count += count_neighbors(neighbor_OH, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_asp ) {
					std::string neighbor_OD1 ("OD1") ;
					neighbor_count += count_neighbors(neighbor_OD1, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
					std::string neighbor_OD2 ("OD2") ;
					neighbor_count += count_neighbors(neighbor_OD2, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_ser ) {
					std::string neighbor_OG ("OG") ;
					neighbor_count += count_neighbors(neighbor_OG, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_thr ) {
					std::string neighbor_OG1 ("OG1") ;
					neighbor_count += count_neighbors(neighbor_OG1, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_glu ) {
					std::string neighbor_OE1 ("OE1") ;
					neighbor_count += count_neighbors(neighbor_OE1, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
					std::string neighbor_OE2 ("OE2") ;
					neighbor_count += count_neighbors(neighbor_OE2, pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
				if ( pose.residue(res_count_neighbor).aa() == core::chemical::AA::aa_gln ) {
					std::string neighbor_OE1 ("OE1") ;
					neighbor_count += count_neighbors(neighbor_OE1, pose,distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness, res_count_target, target_atom_vector_start, distance_internal, target_vector, res_count_neighbor);
				}
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}


// Cone FA Closest Atom Method
void cone_FA_nc_anyatom_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness, core::Real angle_fn_midpoint, core::Real angle_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		std::string target_atom_vector_start ("CA");
		std::string target_atom ("CB");
		if ( pose.residue(res_count_target).aa() == core::chemical::AA::aa_gly ) {
			target_atom = "1HA";
		}
		numeric::xyzVector<core::Real> target_vector = pose.residue(res_count_target).xyz(target_atom) - pose.residue(res_count_target).xyz(target_atom_vector_start);
		core::Real const distance_internal = pose.residue(res_count_target).xyz(target_atom).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				std::string neighbor_atom (pose.residue(res_count_neighbor).atom_name(1));
				for ( core::Size atom_count = 2; atom_count <= pose.residue(res_count_neighbor).natoms(); atom_count++ ) {
					if ( pose.residue(res_count_neighbor).xyz(atom_count).distance(pose.residue(res_count_target).xyz(target_atom)) < pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom)) ) {
						neighbor_atom = pose.residue(res_count_neighbor).atom_name(atom_count);
					}
				}
				numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz(neighbor_atom) - pose.residue(res_count_target).xyz(target_atom_vector_start);
				core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(pose.residue(res_count_target).xyz(target_atom_vector_start));
				numeric::xyzVector<core::Real> norm_target_vector = target_vector/distance_internal;
				numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
				core::Real angle = std::acos(norm_target_vector.dot(norm_neighbor_vector));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}
// Cone Centroid Method
void cone_cen_nc_calculator ( pose::Pose pose, core::Real distance_fn_midpoint, core::Real distance_fn_steepness, core::Real angle_fn_midpoint, core::Real angle_fn_steepness ) {
	std::ostringstream output;
	output << "Residue_Type\tResidue_Number\tNeighbor_Count" << std::endl;
	core::Size const num_residues (pose.total_residue());
	std::string target_vector_start ("CA");
	std::string target ("CEN");
	for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
		core::Real neighbor_count (0.0);
		numeric::xyzVector<core::Real> target_vector = pose.residue(res_count_target).xyz(target) - pose.residue(res_count_target).xyz(target_vector_start);
		core::Real const distance_internal = pose.residue(res_count_target).xyz(target).distance(pose.residue(res_count_target).xyz(target_vector_start));
		for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
			if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
				numeric::xyzVector<core::Real> neighbor_vector = pose.residue(res_count_neighbor).xyz("CEN") - pose.residue(res_count_target).xyz(target_vector_start);
				core::Real const distance_to_neighbor = pose.residue(res_count_neighbor).xyz("CEN").distance(pose.residue(res_count_target).xyz(target_vector_start));
				numeric::xyzVector<core::Real> norm_target_vector = target_vector/distance_internal;
				numeric::xyzVector<core::Real> norm_neighbor_vector = neighbor_vector/distance_to_neighbor;
				core::Real angle = std::acos(norm_target_vector.dot(norm_neighbor_vector));
				neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance_to_neighbor-distance_fn_midpoint)))*1.0/(1.0 + std::exp(angle_fn_steepness*(angle-angle_fn_midpoint)));
			}
		}
		output << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
	}
	output_results(output.str());
}



// Main program
int main( int argc, char * argv [] ){
	try {

		NEW_OPT( dist_midpoint, "midpoint of distance falloff", 9.0);
		NEW_OPT( dist_steepness, "exponent factor for distance falloff", 1.0);
		NEW_OPT( angle_midpoint, "midpoint of angle falloff", numeric::NumericTraits< float >::pi()/2.0);
		NEW_OPT( angle_steepness, "exponent factor for angle falloff", numeric::NumericTraits< float >::pi()*2.0);
		const std::string neighbor_method_cone = "cone";
		const std::string neighbor_method_sphere = "sphere";
		OPT( solvent_exposure::method );
		NEW_OPT( centroid_version, "centroid versions of the neighbor counts", false);
		NEW_OPT( neighbor_closest_atom, "for FA neighbor count, use closest atom per residue as neighbor", false);
		NEW_OPT( bbamide, "For using the NH Vector instead of the CA-CB vector", false);


		// Initialize Rosetta
		devel::init( argc, argv );

		// Namespaces
		using namespace core;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace pose;

		core::Real distance_fn_midpoint = option[ dist_midpoint ];
		core::Real distance_fn_steepness = option[ dist_steepness ];
		core::Real angle_fn_midpoint = option[ angle_midpoint ];
		core::Real angle_fn_steepness = option[ angle_steepness ];
		std::string neighbor_method = option[ solvent_exposure::method ].value();
		bool centroid_pose = option[ centroid_version ];
		bool closest_atom = option[ neighbor_closest_atom ];
		bool amide = option[ bbamide ];
		if ( !option[in::file::s].user() ) {
			utility_exit_with_message("You specify a pdb file using the flag '-in:file:s'!");
		}
		std::string pdb_file = option[in::file::s]()[1];
		pose::Pose pose;
		core::import_pose::pose_from_file(pose,pdb_file);

		// Calculate NC's for Centroid
		if ( centroid_pose == true ) {
			if ( !option[in::file::centroid].user() ) {
				utility_exit_with_message("You need to use the flag '-in:file:centroid' if you wish to find the centroid based solvent exposure");
			}
			if ( neighbor_method == neighbor_method_sphere ) {
				sphere_cen_nc_calculator( pose, distance_fn_midpoint, distance_fn_steepness );
			} else /*if ( neighbor_method_vector[1] == neighbor_method_cone )*/ {
				cone_cen_nc_calculator( pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness );
			}
		} else {
			// Calculate NC's for FA (centroid_pose == false)
			if ( neighbor_method == neighbor_method_sphere ) {
				if ( closest_atom == true ) {
					sphere_FA_nc_anyatom_calculator( pose, distance_fn_midpoint, distance_fn_steepness );
				} else /*if ( closest_atom == false )*/ {
					sphere_FA_nc_calculator( pose, distance_fn_midpoint, distance_fn_steepness );
				}
			} else /*if ( neighbor_method_vector[1] == neighbor_method_cone )*/ {
				if ( closest_atom == true ) {
					cone_FA_nc_anyatom_calculator( pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness );
				} else /*if ( closest_atom == false )*/ {
					if ( amide == true ) {
						cone_FA_nc_calculator_amide( pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness );
					} else {
						cone_FA_nc_calculator( pose, distance_fn_midpoint, distance_fn_steepness, angle_fn_midpoint, angle_fn_steepness );
					}
				}
			}
		}
	}



catch ( utility::excn::Exception const & e ) {
	e.display();
	return -1;
}

	return 0;

} //main
