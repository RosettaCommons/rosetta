// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief declaration generic functions returning info on a constellation
/// @author Andrea Bazzoli

#ifndef INCLUDED_cnl_info_hh
#define INCLUDED_cnl_info_hh

#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

using core::pose::Pose;
using numeric::xyzVector;
using utility::vector1;
using core::Size;
using core::Real;

namespace devel {
namespace constel {

/// @brief returns the center of mass of a constellation
xyzVector<Real> cnl_com(vector1<Size> const &cnl, Pose const &ps);

} // constel
} // devel

#endif

