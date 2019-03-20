// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/hashing/MinimalClashHash.fwd.hh
/// @brief  Minimal hashing clash to find clashes as fast as possible
/// @details A rather limiting design choice of this class is to assume
///          that all atoms have the same radius. xyzStripeHash is
///          the accurate version of this class.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_numeric_VoxelArray_fwd_hh
#define INCLUDED_numeric_VoxelArray_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <numeric/types.hh>

namespace numeric {


template< class _Float=float, class _Value=float >
class VoxelArray;

// I think we're better off not making templated typedefs...
// typedef utility::pointer::shared_ptr< VoxelArray > VoxelArrayOP;
// typedef utility::pointer::shared_ptr< VoxelArray const > VoxelArrayCOP;
// typedef utility::pointer::weak_ptr< VoxelArray const > VoxelArrayCAP;

}

#endif
