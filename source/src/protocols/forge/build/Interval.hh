// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/build/instructions/Interval.fwd.hh
/// @brief  simple struct defining a closed interval of residues
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_Interval_hh
#define INCLUDED_protocols_forge_build_Interval_hh

// unit headers
#include <protocols/forge/build/Interval.fwd.hh>

// project headers
#include <core/types.hh>

// C++ headers
#include <utility/assert.hh>

namespace protocols {
namespace forge {
namespace build {


/// @brief simple struct defining a closed interval of residues [left, right]
///  where left <= right
struct Interval {


	typedef core::Size Size;


	/// @brief default constructor
	inline
	Interval() :
		left( 0 ),
		right( 0 )
	{}


	/// @brief value constructor
	inline
	Interval(
		Size const l,
		Size const r
	) :
		left( l ),
		right( r )
	{
		assert( left <= right );
	}


	/// @brief copy constructor
	inline
	Interval( Interval const & rval ) :
		left( rval.left ),
		right( rval.right )
	{}


	/// @brief default destructor
	inline
	~Interval() {}


	/// @brief copy assignment
	inline
	Interval & operator =( Interval const & rval ) {
		if ( this != &rval ) {
			left = rval.left;
			right = rval.right;
		}
		return *this;
	}


	/// @brief operator <, lexicographic ordering
	inline
	bool operator <( Interval const & rval ) const {
		return (
			( left < rval.left ? true :
			( rval.left < left ? false : // left == rval.left
			( right < rval.right ) ) )
		);
	}


	/// @brief operator ==
	inline
	bool operator ==( Interval const & rval ) const {
		return ( left == rval.left && right == rval.right );
	}


	/// @brief length of the interval
	inline
	Size length() const {
		return right - left + 1;
	}


	/// @brief do the two intervals intersect?
	inline
	bool intersects( Interval const & rval ) const {
		return !( left > rval.right || rval.left > right );
	}


	/// @brief is a point contained in this interval?
	inline
	bool contains( Size const point ) const {
		return ( left <= point && point <= right );
	}


	/// @brief left endpoint
	Size left;


	/// @brief right endpoint
	Size right;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_instructions_Interval_HH */
