// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util/burial_utilities.cc
/// @brief Utilities for determining whether a point is within a cone.
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals; added determine_whether_point_is_buried() function.)

#include <core/select/util/burial_utilities.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

static basic::Tracer TR( "core.select.util.burial_utilities" );

namespace core {
namespace select {
namespace util {

#define NFOLD_REDUCTION 10 //A tenfold reduction in the distance function value is the point at which we no longer consider whether we might be in the distance cone.


utility::vector1< core::Real >
calc_sc_neighbors( pose::Pose const & pose,
	core::Real const angle_exponent /* 2.0 */,
	core::Real const angle_shift_factor /* 0.5 */,
	core::Real const dist_exponent /* 1.0 */,
	core::Real const dist_midpoint /* 9.0 */,
	core::Real const rsd_neighbor_denominator /* 1.0 */,
	bool const asu_only /* false */ )
{
	utility::vector1< Real > rsd_sc_neighbors;

	core::Size target_pose_size;

	//check for symmetry
	bool sym_check = core::pose::symmetry::is_symmetric( pose );
	if ( sym_check && asu_only ) {
		TR.Info << "Pose is symmetric and asu_only=True. Calculating sc_neighbors for asymmetric unit residues ONLY, but it is calculated against the symmetric pose." << std::endl;
		target_pose_size = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		TR.Debug << "Pose.size() = " << pose.size() << std::endl;
		TR.Debug << "ASU size = " << target_pose_size << std::endl;
	} else {
		target_pose_size = pose.size();
	}

	for ( Size i = 1; i <= target_pose_size; ++i ) {
		Real my_neighbors(0.0);

		numeric::xyzVector< Real > my_sc_coordinates;
		numeric::xyzVector< Real > my_bb_coordinates;

		if ( pose.residue( i ).name3() == "GLY" ) {
			my_sc_coordinates = pose.residue(i).atom(pose.residue(i).atom_index("2HA")).xyz() ;
			my_bb_coordinates = pose.residue(i).atom(pose.residue(i).atom_index("CA")).xyz() ;
		} else {
			if ( pose.residue(i).is_polymer() ) {
				my_sc_coordinates = pose.residue(i).atom(pose.residue(i).first_sidechain_atom()).xyz() ;
				core::Size parent_atom_index = pose.residue(i).icoor( pose.residue(i).first_sidechain_atom() ).stub_atom1().atomno();
				my_bb_coordinates = pose.residue(i).atom( parent_atom_index ).xyz() ;
			} else {
				rsd_sc_neighbors.push_back(0); //For now, ligands do not have their neighbours counted.  This could change in the future.
				continue;
			}
		}

		numeric::xyzVector< Real > my_sc_vector = (my_sc_coordinates - my_bb_coordinates).normalize() ;

		for ( Size j = 1; j <= pose.size(); ++j ) {

			//TR.Debug << "Calculating resi " << i << " against resi " << j << std::endl;

			if ( i != j ) {

				numeric::xyzVector< Real > other_bb_coordinates;
				if ( pose.residue(j).name3() == "GLY" ) {
					other_bb_coordinates = pose.residue(j).atom(pose.residue(j).atom_index("CA")).xyz();
				} else {
					if ( pose.residue(j).is_polymer() ) { //If this is a polymer atom, use the parent of the first sidechain atom.
						core::Size parent_atom_index = pose.residue(j).icoor( pose.residue(j).first_sidechain_atom() ).stub_atom1().atomno();
						other_bb_coordinates = pose.residue(j).atom( parent_atom_index ).xyz();
					} else { //If this is not a polymer residue, use the nbr_atom:
						core::Size nbr_atom_index = pose.residue(j).nbr_atom();
						other_bb_coordinates = pose.residue(j).atom( nbr_atom_index ).xyz();
					}
				}

				my_neighbors += calculate_point_in_cone( other_bb_coordinates, my_sc_vector, my_sc_coordinates, angle_exponent, angle_shift_factor, dist_exponent, dist_midpoint );
			}
		}
		rsd_sc_neighbors.push_back(my_neighbors / rsd_neighbor_denominator);
		TR.Debug << "Resi " << i << "; my_neighbors: " << my_neighbors << std::endl;
	}

	return rsd_sc_neighbors;
}

/// @brief Given a point in 3D space, determine whether or not that point is buried by the method of sidechain neighbor cones.
/// @note A crude distance cutoff metric is also used to make this calculation more efficient.
/// @details Returns true for burial, false otherwise.  Updated 1 March 2019 to ignore ligands.
/// @param[in] point_coordinates The 3D coordinates of the point in space.
/// @param[in] pose The pose, for reference.
/// @param[in] angle_exponent A value defining how rapidly the cone falls off in the angular direction.
/// @param[in] angle_shift_factor A value that shifts the angluar falloff.
/// @param[in] dist_exponent A value defining how rapidly the cone falls off with distance.
/// @param[in] dist_midpoint A value defining the midpoint of the distance falloff.
/// @param[in] burial_threshold The cutoff for considering a point to be buried.  Roughly, this is the number of cones that the point must be inside
/// in order to consider it "buried".
bool determine_whether_point_is_buried (
	numeric::xyzVector<core::Real> const &point_coordinates,
	core::pose::Pose const &pose,
	core::Real const angle_exponent,
	core::Real const angle_shift_factor,
	core::Real const dist_exponent,
	core::Real const dist_midpoint,
	core::Real const burial_threshold
) {
	core::Real accumulator(0.0);
	core::Real const dist_cutoff_sq( std::pow(std::log( NFOLD_REDUCTION - 1 ) + dist_exponent*dist_midpoint, 2) ); //Used to skip calculation for distant points.

	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		core::chemical::ResidueType const &restype( pose.residue_type(ir) );
		if ( restype.is_ligand() ) continue;
		core::id::AtomID at2( restype.aa() == core::chemical::aa_gly ? restype.atom_index("2HA") : restype.first_sidechain_atom(), ir );
		core::id::AtomID at1( restype.icoor(at2.atomno()).stub_atom1().atomno(), ir );
		numeric::xyzVector< core::Real > const conevect_coordinate_2( pose.xyz(at2) );

		//Some crude distance checks to accelerate this:
		if ( point_coordinates.distance_squared( conevect_coordinate_2 ) > dist_cutoff_sq ) continue;

		//Do the actual calculation:
		numeric::xyzVector< core::Real > const conevect( conevect_coordinate_2 - pose.xyz(at1) );
		accumulator += calculate_point_in_cone( point_coordinates, conevect, conevect_coordinate_2, angle_exponent, angle_shift_factor, dist_exponent, dist_midpoint );
		if ( accumulator > burial_threshold ) return true;
	}
	return false;
}

} //namespace util
} //namespace select
} //namespace core
