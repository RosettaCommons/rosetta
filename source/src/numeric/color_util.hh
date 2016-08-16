// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/color_util.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_color_util_HH
#define INCLUDED_numeric_color_util_HH

#include <platform/types.hh>
#include <numeric/xyzVector.fwd.hh>
namespace numeric {

/// @brief convert an RGB color to HSV
numeric::xyzVector<platform::Real> rgb_to_hsv(platform::Real r,platform::Real b, platform::Real g);

/// @brief convert and RGB color to HSV
numeric::xyzVector<platform::Real> rgb_to_hsv(numeric::xyzVector<platform::Real> rgb_triplet);

/// @brief convert an HSV color to RGB
numeric::xyzVector<platform::Real> hsv_to_rgb(platform::Real h, platform::Real s, platform::Real v);

/// @brief convert an HSV color to RGB
numeric::xyzVector<platform::Real> hsv_to_rgb(numeric::xyzVector<platform::Real> hsv_triplet);


}

#endif /* COLOR_UTIL_HH */
