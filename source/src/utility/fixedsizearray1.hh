// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief unresizable vector whose size is known at compile time,
/// which may be allocated on the stack, and which indexes from 1.

#ifndef INCLUDED_utility_fixedsizearray1_hh
#define INCLUDED_utility_fixedsizearray1_hh

// Unit headers
#include <utility/fixedsizearray1.fwd.hh>

// C++ headers
#include <array>
#include <utility/assert.hh>

namespace utility {

/// Requirements:
/// S must be a positive integer
/// T must be default constructable

template < typename T, platform::Size S >
class fixedsizearray1
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
	/// Constructors

	fixedsizearray1() {
		array_.fill( value_type{} ); // Zero initialize for POD types, call the default constructor for classes.
	}

	fixedsizearray1( value_type const & def ) {
		array_.fill( def );
	}

	/// @brief initializer list construction
	fixedsizearray1( std::initializer_list< T > const & source ) {
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
	fixedsizearray1< T, S > const &
	operator = ( value_type const & val ) {
		array_.fill( val );
		return *this;
	}

public:
	/// Mutators and accessors

	value_type &
	operator [] ( Size index ) {
		assert( range( index ) ); // debug_assert() gives compile errors for gcc 4.8 release_debug compile
		int zero_based_index( index - 1 ); // Needed for gcc 4.9 static compilation, as Size as an index throws a warning
		return array_[ zero_based_index ];
	}

	value_type const &
	operator [] ( Size index ) const  {
		assert( range( index ) ); // debug_assert() gives compile errors for gcc 4.8 release_debug compile
		int zero_based_index( index - 1 ); // Needed for gcc 4.9 static compilation, as Size as an index throws a warning
		return array_[ zero_based_index ];
	}

	bool
	operator == ( fixedsizearray1< T, S > const & other ) const {
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
		return index > (Size) 0 && index <= (Size) S;
	}

private:
	/// Data
	std::array< T, S > array_;
};


} // namespace utility

#endif
