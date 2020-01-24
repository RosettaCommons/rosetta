// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief unresizable vector whose size is known at compile time,
/// which may be allocated on the stack, and which indexes from 0.

#ifndef INCLUDED_utility_fixedsizearray0_hh
#define INCLUDED_utility_fixedsizearray0_hh

// Unit headers
#include <utility/fixedsizearray0.fwd.hh>

// C++ headers
#include <array>
#include <utility/assert.hh>

namespace utility {

/// Requirements:
/// S must be a positive integer
/// T must be default constructable.

template < typename T, platform::Size S >
class fixedsizearray0
{
public:
	typedef platform::Size Size;
	typedef platform::SSize SSize;

	typedef T value_type;
	typedef typename std::array<T,S>::iterator iterator;
	typedef typename std::array<T,S>::const_iterator const_iterator;

	//typedef  typename T &        reference;
	//typedef  typename T const &  const_reference;

public:
	/// Constructors and the assigmnet operator

	fixedsizearray0() {
		array_.fill( value_type{} ); // Zero initialize for POD types, call the default constructor for classes.
	}

	fixedsizearray0( value_type const & def ) {
		array_.fill( def );
	}

	/// @brief initializer list construction
	fixedsizearray0( std::initializer_list< T > const & source ) {
		if ( source.size() == 1 ) {
			// If we just pass one element, fill the array with it.
			array_.fill( *source.begin() );
		} else {
			assert( S == source.size() );
			Size pos = 0;
			for ( auto const & val: source ) {
				array_[pos] = val;
				++pos;
			}
		}
	}

	/// @brief Full-array assignment.
	fixedsizearray0< T, S > const &
	operator = ( value_type const & val ) {
		array_.fill( val );
		return *this;
	}

public:
	/// Mutators and accessors

	value_type &
	operator [] ( Size index ) {
		debug_assert( range( index ) ); // debug_assert() gives compile errors for gcc 4.8 release_debug compile
		return array_[ index ];
	}

	value_type const &
	operator [] ( Size index ) const  {
		debug_assert( range( index ) ); // debug_assert() gives compile errors for gcc 4.8 release_debug compile
		return array_[ index ];
	}

	bool
	operator == ( fixedsizearray0< T, S > const & other ) const {
		return array_ == other.array_;
	}


	Size
	size() const {
		return S;
	}

public:
	/// Iterators
	iterator
	begin() {
		return array_.begin();
	}

	iterator
	end() {
		return array_.end();
	}

	const_iterator
	begin() const {
		return array_.cbegin();
	}

	const_iterator
	end() const {
		return array_.cend();
	}

protected:

	bool
	range( Size index ) const {
		//return index > (Size) 0 && index <= (Size) S;
		return index >= (Size) 0 && index < (Size) S; // should be indexed by 0, changed by Georg Kuenze 05/10/2016
	}

private:
	/// Data
	std::array< T, S > array_;
};


} // namespace utility

#endif
