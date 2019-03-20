// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/hashing/MinimalClashHash.fwd.hh
/// @brief  Minimal hashing class to find clashes as fast as possible
/// @details A rather limiting design choice of this class is to assume
///          that all atoms have the same radius. xyzStripeHash is
///          the accurate version of this class.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_numeric_geometry_hashing_MinimalClashHash_fwd_hh
#define INCLUDED_numeric_geometry_hashing_MinimalClashHash_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <numeric/types.hh>

namespace numeric {
namespace geometry {
namespace hashing {

class MinimalClashHash;
typedef utility::pointer::shared_ptr< MinimalClashHash > MinimalClashHashOP;
typedef utility::pointer::shared_ptr< MinimalClashHash const > MinimalClashHashCOP;
typedef utility::pointer::weak_ptr< MinimalClashHash const > MinimalClashHashCAP;

}
}
}

#endif
