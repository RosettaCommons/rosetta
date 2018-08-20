// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/atomic_depth/util.hh
/// @brief  Util functions for atomic_depth calculations.
/// @author Brian Coventry

#ifndef INCLUDED_core_scoring_atomic_depth_util_hh
#define INCLUDED_core_scoring_atomic_depth_util_hh

#include <core/scoring/atomic_depth/AtomicDepth.hh>

#include <core/id/AtomID_Map.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace atomic_depth {


/// @brief Calculate depth of all atoms from edge of Sasa surface.
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	Real probe_radius = 1.4,
	bool poly_leu_depth = false
);

/// @brief Calculate depth of all atoms from edge of Sasa surface.
/// @detail If atomic_depth is a nullptr, will create one and return it.
/// @detail  Note that probe_radius and poly_leu_depth are ignored if atomic_depth is not null
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	AtomicDepthOP & depth,
	Real probe_radius = 1.4,
	bool poly_leu_depth = false
);

/// @brief Calculate depth of all atoms from edge of Sasa surface.
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	core::id::AtomID_Map< bool > depth_atoms,
	Real probe_radius = 1.4,
	bool poly_leu_depth = false
);

/// @brief Calculate depth of all atoms from edge of Sasa surface.
/// @detail If atomic_depth is a nullptr, will create one and return it.
/// @detail  Note that probe_radius and poly_leu_depth are ignored if atomic_depth is not null
core::id::AtomID_Map< core::Real >
atomic_depth(
	pose::Pose const & pose,
	core::id::AtomID_Map< bool > depth_atoms,
	AtomicDepthOP & depth,
	Real probe_radius = 1.4,
	bool poly_leu_depth = false
);

/// @brief Find all atoms deeper than threshold from the edge of the sasa surface.
core::id::AtomID_Map< bool >
atoms_deeper_than(
	pose::Pose const & pose,
	Real threshold,
	bool invert = false,
	Real probe_radius = 1.4,
	bool poly_leu_depth = false
);




} // namespace atomic_depth
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_atomic_depth_AtomicDepth_fwd_hh
