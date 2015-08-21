// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_utility_FixedSizeLexicographicalIterator_tmpl_hh
#define INCLUDED_utility_FixedSizeLexicographicalIterator_tmpl_hh

// Unit headers
#include <utility/FixedSizeLexicographicalIterator.hh>

// C++ headers
///#include <iostream> // debug only

namespace utility {

/// @details -- fake that there's one alphabet and that it has size 0.
template < platform::Size S >
FixedSizeLexicographicalIterator< S >::FixedSizeLexicographicalIterator() :
dim_sizes_( 0 ),
curr_pos_( 1 )
{
}

/// @details No alphabet may have size 0.
template < platform::Size S >
FixedSizeLexicographicalIterator< S >::FixedSizeLexicographicalIterator(
	fsarray const & dim_sizes
) :
dim_sizes_( dim_sizes ),
curr_pos_( 1 )
{
	for ( Size ii = 1; ii <= S; ++ii ) {
		assert ( dim_sizes[ ii ] > 0 );
	}
}

template < platform::Size S >
void FixedSizeLexicographicalIterator< S >::set_dimension_sizes(
	fsarray const & dim_sizes
)
{
	dim_sizes_ = dim_sizes;
	for ( Size ii = 1; ii <= S; ++ii ) {
		debug_assert( dim_sizes[ ii ] > 0 );
	}
	begin();
}


/// @details First element is the 1 string
template < platform::Size S >
void
FixedSizeLexicographicalIterator< S >::begin()
{
	std::fill( curr_pos_.begin(), curr_pos_.end(), 1 );
}

/// @details both curr_pos_ and dim_sizes_ must have size greater than 0.
/// Therefore, even the empty set of alphabets (used by the default ctor)
/// must contain 1 alphabet.
template < platform::Size S >
bool FixedSizeLexicographicalIterator< S >::at_end() const
{
	//for ( Size ii = 1; ii <= dim_sizes_.size(); ++ii ) {
	// if ( curr_pos_[ ii ] != dim_sizes_[ ii ] ) return false;
	//}
	//return true;
	return curr_pos_[ 1 ] > dim_sizes_[ 1 ];
}

/// @details Return the number of dimensions advanced by the increment operator
/// so that an outside observer can track the progress of the iterator
template < platform::Size S >
typename FixedSizeLexicographicalIterator< S >::Size
FixedSizeLexicographicalIterator< S >::operator ++ () {
	// do not increment past the end position, return 0 to say no dimensions were advanced.
	if ( at_end() ) return 0;

	for ( Size ii = dim_sizes_.size(), n_advanced = 1; ii > 1; --ii, ++n_advanced ) {
		if ( curr_pos_[ ii ] == dim_sizes_[ ii ] ) {
			curr_pos_[ ii ] = 1;
		} else {
			++curr_pos_[ ii ];
			return n_advanced;
		}
	}
	++curr_pos_[ 1 ];
	return S; // all dimensions were advanced
}

/// @brief Give an integer index of the current state.  This can be used to reset the
/// lexicographical iterator to the current state again later.
template < platform::Size S >
typename FixedSizeLexicographicalIterator< S >::Size
FixedSizeLexicographicalIterator< S >::index() const {
	Size ind( 1 );  /// index from 1 to numstates, instead of from 0 to numstates - 1;
	Size dimprods( 1 );
	for ( Size ii = S; ii >= 1; --ii ) {
		ind += ( curr_pos_[ ii ] - 1 ) * dimprods;
		dimprods *= dim_sizes_[ ii ];
	}
	return ind;
}


/// @details No requirement that the index doesn't "overflow" the iterator; however,
/// if it should overflow, then at_end() will return true.
template < platform::Size S >
void
FixedSizeLexicographicalIterator< S >::set_position_from_index( Size index ) {
	index -= 1;

	Size dimprods( 1 );
	for ( Size ii = 2; ii <= S; ++ii ) {
		dimprods *= dim_sizes_[ ii ];
	}
	///std::cout << "start: dimprods " << dimprods << " index: " << index << std::endl;
	/// loop over all but the last dimension;
	for ( Size ii = 1; ii < S; ++ii ) {
		curr_pos_[ ii ] = 1 + ( index / dimprods ); // integer division
		index = index % dimprods;
		dimprods /= dim_sizes_[ ii + 1 ];

		//std::cout << ii << ": dimprods " << dimprods << " curr_pos_[ ii ] " << curr_pos_[ ii ] << " index: " << index << std::endl;

	}
	curr_pos_[ dim_sizes_.size() ] = 1 + index;
}

/// @details set the higher dimensions to their maximum values, and simply invoke
/// the ++operator which will handle the logic of rolling over lower dimensions
/// if the nth dimension is at its maximum value.
template < platform::Size S >
typename FixedSizeLexicographicalIterator< S >::Size
FixedSizeLexicographicalIterator< S >::continue_at_dimension( Size dim ) {
	for ( Size ii = dim+1; ii <= dim_sizes_.size(); ++ii ) {
		curr_pos_[ ii ] = dim_sizes_[ ii ];
	}
	return ++(*this);
}

template < platform::Size S >
typename FixedSizeLexicographicalIterator< S >::Size
FixedSizeLexicographicalIterator< S >::num_states_total() const {
	Size total = 1;
	for ( Size ii = 1; ii <= dim_sizes_.size(); ++ii ) {
		total *= dim_sizes_[ ii ];
	}
	return total;
}


}

#endif
