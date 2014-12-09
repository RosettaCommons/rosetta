// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/FixedSizeLexicographicalIterator.hh
/// @brief  Class for iterating across all words of a compile-time-known fixed-length, composed of letters from
/// varying alphabets, where each alphabet is represented by its size.  This iteration is
/// performed in lexicographical order.  All allocation of a FixedSizeLexicographicalIterator takes
/// place on the stack, and is therefore much faster than a regular LexigographicalIterator.
///
/// @author  Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_FixedSizeLexicographicalIterator_hh
#define INCLUDED_utility_FixedSizeLexicographicalIterator_hh

// Unit headers
#include <utility/FixedSizeLexicographicalIterator.fwd.hh>
#include <utility/fixedsizearray1.hh>

namespace utility {

template < platform::Size S >
class FixedSizeLexicographicalIterator {
public:
	typedef platform::Size Size;
	typedef fixedsizearray1< Size, S > fsarray;

public:
	/// @brief default constructor -- requires that set_dimension_sizes later be called.
	/// or the default alphabet (size 0).
	FixedSizeLexicographicalIterator();

	/// @brief Constructor with dimension sizes.
	FixedSizeLexicographicalIterator( fsarray const & dim_sizes );

	/// @brief Set the dimension sizes; resets the curr_pos_ to the beginning of the
	/// enumeration.
	void set_dimension_sizes( fsarray const & dim_sizes );

	/// @brief reset the iterator to the beginning string (1,1,1,...)
	void begin();

	/// @brief Is the iterator at the end?
	bool at_end() const;

	/// @brief Increment the iterator and return the number of dimensions that were advanced.
	/// the number of advanced dimensions ranges from 0 to ndims.  0 is returned only if
	/// the iterator is at the end.
	Size operator ++ ();

	/// @brief The number of dimensions
	Size size() const {
		return dim_sizes_.size();
	}

	/// @brief The number of dimensions
	Size ndims() const {
		return size();
	}

	/// @brief Access the ith dimension (from most-significant digit to least).
	/// If the iterator pointed to the string "MACE", then dimension "2" refers to the
	/// position holding "A". Unsigned dimension input.
	inline
	Size operator[] ( Size dim ) const {
		return curr_pos_[ dim ];
	}

	/// @brief
	inline
	Size dimsize( Size dim ) const {
		return dim_sizes_[ dim ];
	}


	/// @brief Give an integer index of the current state.  This can be used to reset the
	/// lexicographical iterator to the current state again later.
	Size index() const;

	/// @brief Set the state of the lexicographical iterator using a particular index.
	void set_position_from_index( Size index );

	/// @brief Advance the nth dimension to its next value and reset the higher dimensions
	/// to their initial values.  E.g. If there were four dimensions of size 5, and
	/// the current state was [ 1, 3, 2, 4 ], then continue_at_dimension( 2 ) would result
	/// in the state [ 1, 4, 1, 1 ], and if the state were   [ 1, 5, 2, 3 ], then
	/// continue_at_dimension( 2 ) would result in the state [ 2, 1, 1, 1 ].
	/// Returns the number of dimensions that were advanced (0 if at_end, but otherwise, >= dim)
	Size continue_at_dimension( Size dim );

	/// @brief Returns the number of states that could be enumerated
	Size num_states_total() const;

private:
	fsarray dim_sizes_;
	fsarray curr_pos_;

};

}

#endif
