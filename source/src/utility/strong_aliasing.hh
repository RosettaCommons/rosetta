// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/strong_aliasing.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_utility_strong_aliasing_hh
#define INCLUDED_utility_strong_aliasing_hh

#include <platform/types.hh>

namespace utility {

//#include <NamedType/named_type.hpp>

/*
example usage:

OLD:
void draw_rectangle( int width, int height );

void foo( int h, int w ){
draw_rectangle( h, w ); //logic error, may not be caught at runtime
draw_rectangle( w, h ); //correct
}


NEW:
using Width = utility::StrongSize< struct Width_ >;
using Height = utility::StrongSize< struct Height_ >;

void draw_rectangle( Width width, Height height );

void foo( Size h, Size w ){
draw_rectangle( h, w ); //compiler error
draw_rectangle( w, h ); //compiler error
draw_rectangle( Width( h ), Height( w ) ); //logic error, but the developer should be able to notice that this is wrong
draw_rectangle( Width( w ), Height( h ) ); //correct
}

*/


//This might have some overhead for things like strings and vectors that have heap data.
template< typename T, typename Key >
struct StrongT {
	T value;

	explicit StrongT( T v ):
		value( v )
	{}

	operator T() const {
		return value;
	}

};

template< typename Key >
using StrongReal = StrongT< platform::Real, Key >;

template< typename Key >
using StrongSize = StrongT< platform::Size, Key >;

}  // namespace utility

#endif  // INCLUDED_utility_strong_aliasing_HH
