// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/Hit.fwd.hh
/// @brief  Hit class forward declaration
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_Hit_fwd_hh
#define INCLUDED_protocols_match_Hit_fwd_hh

// Project headers
#include <core/types.hh>

// Numeric headers
#include <numeric/kdtree/WrappedPrimitive.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <list>

#include <utility/fixedsizearray1.fwd.hh>


namespace protocols {
namespace match {

typedef numeric::geometry::hashing::Size3 Size3;
typedef numeric::geometry::hashing::Size4 Size4;
typedef numeric::geometry::hashing::Real6 Real6;

/// A hit is defined as the integer representation of the upstream partner,
/// the integer representation of the downstream partner's internal geometry,
/// and the 6 DOFs describing the placement of the downstream partner in
/// the global coordinate frame.
///
/// Hit hit;
/// hit.first[ 1 ] == scaffold build id
/// hit.first[ 2 ] == upstream conf id
/// hit.first[ 3 ] == which external geometry generated this hit
/// hit.first[ 4 ] == downstream conf id
/// hit.second     == 6 DOFs describing the downstream conf's rigid body conformation.

///typedef std::pair< Size4, Real6 > Hit;

class Hit;

typedef utility::vector1< Hit > match;

typedef numeric::kdtree::WrappedPrimitive< std::list< Hit const * > > HitPtrList;
typedef utility::pointer::shared_ptr< HitPtrList > HitPtrListOP;
typedef utility::pointer::shared_ptr< HitPtrList const > HitPtrListCOP;

class upstream_hit;
class downstream_hit;
struct match_dspos1;

}
}

#endif
