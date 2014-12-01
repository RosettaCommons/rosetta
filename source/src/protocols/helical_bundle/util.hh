// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/util.cc
/// @brief  Headers for utility functions for helical bundle construction.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_util_hh
#define INCLUDED_protocols_helical_bundle_util_hh

// Unit Headers
#include <protocols/moves/Mover.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

	void write_minor_helix_params (
		std::string const &filename,
		utility::vector1 < core::Real > const &r1,
		core::Real const &omega1,
		core::Real const &z1,
		utility::vector1 < core::Real > const &delta_omega1,
		utility::vector1 < core::Real > const &delta_z1
	);

	void read_minor_helix_params (
		std::string const &filename,
		utility::vector1 < core::Real > &r1,
		core::Real &omega1,
		core::Real &z1,
		utility::vector1 < core::Real > &delta_omega1,
		utility::vector1 < core::Real > &delta_z1
	);

	/// @brief Generate the x,y,z coordinates of the mainchain atoms using the Crick equations.
	/// @details Coordinates will be returned as a vector of vectors of xyzVectors.  The outer
	/// index will refer to residue number, and the inner index will refer to atom number.
	/// Returns failed=true if coordinates could not be generated, false otherwise.
	void generate_atom_positions(
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > &outvector,
		core::pose::Pose const &helixpose,
		core::Size const helix_start,
		core::Size const helix_end,
		core::Real const &r0,
		core::Real const &omega0,
		core::Real const &delta_omega0,
		core::Real const &delta_t,
		bool const invert_helix,
		utility::vector1 < core::Real > const &r1,
		core::Real const &omega1,
		core::Real const &z1,
		utility::vector1 < core::Real > const &delta_omega1,
		core::Real const &delta_omega1_all,
		utility::vector1 < core::Real > const &delta_z1,
		bool &failed
	);

	/// @brief Place the helix mainchain atoms based on the Crick equations.
	///
	void place_atom_positions(
		core::pose::Pose &pose,
		utility::vector1 < utility::vector1 < numeric::xyzVector < core::Real >  > > const &atom_positions,
		core::Size const helix_start,
		core::Size const helix_end
	);

	/// @brief Copy backbone bond length values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_bondlengths(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	);

	/// @brief Copy backbone bond angle values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_bondangles(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	);

	/// @brief Copy backbone dihedral values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_dihedrals(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	);

	/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms.
	///
	void align_mainchain_atoms(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	);

	/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms,
	/// moving ONLY the residues involved in the alignment.
	void align_mainchain_atoms_of_residue_range(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	);

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_util_hh
