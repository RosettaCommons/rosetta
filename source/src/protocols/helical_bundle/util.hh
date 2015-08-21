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
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

/// @brief Actual write of the crick_params file data.
/// @details Called by both write_minor_helix_params variants.
/// The outfile ozstream must already be opened, and will not be closed
/// by this function.
void write_crick_params_file_data (
	utility::io::ozstream &outfile,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
);

/// @brief Write out a crick_params file.
/// @details Variant for case of a single-residue repeating unit.
void write_minor_helix_params (
	std::string const &filename,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
);

/// @brief Write out a crick_params file.
/// @details Variant for case of a multi-residue repeating unit.
void write_minor_helix_params (
	std::string const &filename,
	core::Size const &residues_per_repeat,
	utility::vector1 <core::Size> const &atoms_per_residue,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	utility::vector1 < core::Real > const &delta_z1
);

/// @brief Read minor helix parameters from a crick_params file.
///
void read_minor_helix_params (
	std::string const &filename,
	utility::vector1 < core::Real > &r1,
	core::Real &omega1,
	core::Real &z1,
	utility::vector1 < core::Real > &delta_omega1,
	utility::vector1 < core::Real > &delta_z1,
	core::Size &residues_per_repeat,
	utility::vector1 <core::Size> &atoms_per_residue
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
	core::Real const &z1_offset,
	core::Real const &z0_offset,
	bool const invert_helix,
	utility::vector1 < core::Real > const &r1,
	core::Real const &omega1,
	core::Real const &z1,
	utility::vector1 < core::Real > const &delta_omega1,
	core::Real const &delta_omega1_all,
	utility::vector1 < core::Real > const &delta_z1,
	core::Size const residues_per_repeat,
	utility::vector1 <core::Size> const &atoms_per_residue,
	core::Size const repeating_unit_offset, //0 if the first residue is the first residue in the repeating unit; 1 if we're off by 1, etc.
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

/// @brief Given a comma-separated list of residue names, separate these out into a vector of
/// residue names.
/// @details The string_in string is the input; the vect_out vector is the output (which will
/// be reset by this operation).
void parse_resnames( std::string const &string_in, utility::vector1< std::string > &vect_out );

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_util_hh
