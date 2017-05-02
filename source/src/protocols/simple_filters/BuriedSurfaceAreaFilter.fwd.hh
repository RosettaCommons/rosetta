// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedSurfaceAreaFilter.fwd.hh
/// @brief Forward declarations and owning pointer typedefs for the BuriedSurfaceAreaFilter.
/// @details Calculates buried surface area (exposed surface area minus total surface area, on a per-residue basis).  Accepts
/// a residue selector to allow buried subsets to be considered.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_filters {

class BuriedSurfaceAreaFilter; //Forward declaration

//Owning pointers:
typedef utility::pointer::shared_ptr< BuriedSurfaceAreaFilter > BuriedSurfaceAreaFilterOP;
typedef utility::pointer::shared_ptr< BuriedSurfaceAreaFilter const > BuriedSurfaceAreaFilterCOP;

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_fwd_hh
