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

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


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

	StrongT(){}

	explicit StrongT( T v ):
		value( v )
	{}

	operator T() const {
		return value;
	}

	T& operator() () {
		return value;
	}

	T const & operator() () const {
		return value;
	}

	//Operator overloads to make common functions easier:

	StrongT< T, Key > &
	operator+= ( T other ){
		value += other;
		return *this;
	}

	StrongT< T, Key > &
	operator-= ( T other ){
		value -= other;
		return *this;
	}

	//NOTE:
	// + and - operators work implicitly

	//Prefix increment operator
	StrongT< T, Key > &
	operator++ (){
		++value;
		return *this;
	}

	//Prefix decrement operator
	StrongT< T, Key > &
	operator-- (){
		--value;
		return *this;
	}
};



template< typename Key >
using StrongReal = StrongT< platform::Real, Key >;

template< typename Key >
using StrongSize = StrongT< platform::Size, Key >;


//Serialization Macros

#ifdef SERIALIZATION

#define SERIALIZE_STRONG_SIZE_HH( NAME ) using NAME = utility::StrongSize< struct NAME##_ >;\
	template < class Archive > void save( Archive & archive, NAME const & );\
	template < class Archive > void load( Archive & archive, NAME & )

#define SERIALIZE_STRONG_SIZE_CC( NAME ) template < class Archive > void save( Archive & archive, NAME const & t ) { archive( CEREAL_NVP( t.value ) ); } \
	template < class Archive > void load( Archive & archive, NAME & t ) { archive( t.value ); } \
	EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( NAME )


#define SERIALIZE_STRONG_REAL_HH( NAME ) using NAME = utility::StrongReal< struct NAME##_ >;\
	template < class Archive > void save( Archive & archive, NAME const & );\
	template < class Archive > void load( Archive & archive, NAME & )

#define SERIALIZE_STRONG_REAL_CC( NAME ) template < class Archive > void save( Archive & archive, NAME const & t ) { archive( CEREAL_NVP( t.value ) ); } \
	template < class Archive > void load( Archive & archive, NAME & t ) { archive( t.value ); } \
	EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( NAME )

#endif

}  // namespace utility


#endif  // INCLUDED_utility_strong_aliasing_HH
