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

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_util_hh
