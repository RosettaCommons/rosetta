// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/iterators.hh
/// @brief  iterator utilities
/// @details unit tests for these utilities are grouped with core/select/iterators.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_utility_iterators_hh
#define INCLUDED_utility_iterators_hh

// Package headers
#include <platform/types.hh>
#include <utility/vector1.hh>

namespace utility {

class Enumerate1 {
public:
	Enumerate1( platform::Size r ) :
		res_( r )
	{}

	bool operator == ( Enumerate1 const & o ) const { return res_ == o.res_; }
	bool operator != ( Enumerate1 const & o ) const { return res_ != o.res_; }

	platform::Size const & operator * () const { return res_; }
	platform::Size & operator * () { return res_; }

	void operator++(){ ++res_; }
private:
	platform::Size res_ = 0;
};


class SimpleRange1 {
	//see https://godbolt.org/z/YtZVj3
public:
	SimpleRange1( platform::Size n ) :
		nres_( n )
	{}

	Enumerate1 begin() {
		return Enumerate1( 1 );
	}

	Enumerate1 end() {
		return Enumerate1( nres_ + 1 );
	}

private:
	platform::Size nres_ = 0;
};

///@brief iterate from [1, max_inclusive]
///@details USAGE: for( platform::Size const index : enumerate1( vec.size() ) )
inline
SimpleRange1
enumerate1( platform::Size const max_inclusive ){
	return SimpleRange1( max_inclusive );
}

///@brief iterate from [1, c.size()]
///@details USAGE: for( platform::Size const index : indices1( vec1 ) )
template< typename T >
inline
SimpleRange1
indices1( utility::vector1< T > const & v ){
	return enumerate1( v.size() );
}

} // utility namespace

#endif // INCLUDED
