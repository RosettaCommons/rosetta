// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/util.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_util_hh
#define INCLUDED_core_scoring_packstat_util_hh

#include "core/scoring/packstat/types.hh"
#include "core/scoring/packstat/CavityBall.hh"
// AUTO-REMOVED #include "core/scoring/packstat/SimplePDB_Atom.hh"

#include "numeric/xyzVector.hh"


namespace core {
namespace scoring {
namespace packstat {

// searches for sphere/ball in vector sorted by x value
inline std::size_t search_x( Spheres  const & spheres, PackstatReal const x, std::size_t begin, std::size_t end ) {
	if( end - begin < 2 ) return begin;
	size_t mid = ( end - begin ) / 2 + begin;
	if( spheres[mid].xyz.x() <= x ) {
		return search_x( spheres, x, mid,  end  );
	} else {
		return search_x( spheres, x,  1 , mid-1 );
	}
}

inline std::size_t search_x( Spheres  const & spheres, PackstatReal const x ) {
	return search_x( spheres, x, (std::size_t)1, spheres.size() );
}

inline std::size_t search_x( CavBalls const & cbs    , PackstatReal const x, std::size_t begin, std::size_t end ) {
	if( end - begin < 2 ) return begin;
	std::size_t mid = ( end - begin ) / 2 + begin;
	if( cbs[mid].xyz().x() <= x ) {
		return search_x( cbs, x, mid,  end  );
	} else {
		return search_x( cbs, x,  1 , mid-1 );
	}
}

inline std::size_t search_x( CavBalls const & cbs    , PackstatReal const x ) {
	return search_x( cbs, x, (size_t)1, cbs.size() );
}

inline PackstatReal max_rad( Spheres const & s ) {
	PackstatReal maxrad = 0;
	for( SphereCIter i = s.begin(); i != s.end(); ++i ) {
		if( maxrad < i->radius ) maxrad = i->radius;
	}
	return maxrad;
}

inline PackstatReal max_rad( CavBalls const & s ) {
	PackstatReal maxrad = 0;
	for( CavBallCIter i = s.begin(); i != s.end(); ++i ) {
		if( maxrad < i->radius() ) maxrad = i->radius();
	}
	return maxrad;
}

} // namespace packstat
} // namespace scoring
} // namespace core

#endif
