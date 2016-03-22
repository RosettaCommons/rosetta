// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/cryst/CheshireCell.hh
/// @brief A struct used by spacegroup.cc and its helper functions.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_protocols_cryst_CheshireCell_hh
#define INCLUDED_protocols_cryst_CheshireCell_hh

// Rosetta core includes:
#include <core/types.hh>
#include <core/kinematics/RT.fwd.hh>

// Rosetta numeric includes:
#include <numeric/xyzVector.hh>

namespace protocols {
namespace cryst {

/// @brief A struct used by spacegroup.cc and its helper functions.
///
struct CheshireCell {
	CheshireCell() : low(0.0,0.0,0.0),high(0.0,0.0,0.0) {};
	CheshireCell( numeric::xyzVector<core::Real> low_in, numeric::xyzVector<core::Real> high_in ) :
		low( low_in ),
		high( high_in )
	{}
	numeric::xyzVector<core::Real> low, high;
};

}
}

#endif
