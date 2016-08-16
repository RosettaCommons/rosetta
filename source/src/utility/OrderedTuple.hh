// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/OrderedTuple.hh
/// @brief  Class for compairing/sorting data where sort-precidence is
/// in descending order from begin() to end()
/// @author  Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_OrderedTuple_hh
#define INCLUDED_utility_OrderedTuple_hh

// Unit headers
#include <utility/OrderedTuple.fwd.hh>

// Project headers
#include <platform/types.hh>

// C++ headers
#include <utility/assert.hh>

namespace utility {

/// @brief The ordered tuple takes a container class T
/// and defines comparison operators so that the tuple may be sorted.
///
/// @details container class T must define stl-like "const_iterator" and
/// "iterator" typedefs.  T must return iterators with calls to
/// begin() and end(). T must be copyable and assignable.
/// T::value_type must be comparable.
template < class T >
class OrderedTuple {
public:
	typedef platform::Size                     Size;
	typedef T                                  container;
	typedef typename container::const_iterator const_iterator;
	typedef typename container::iterator       iterator;

public:
	/// @brief default constructor
	OrderedTuple() {}

	OrderedTuple( T const & data ) : data_( data ) {}

	void
	assign_data( container const & val ) {
		data_ = val;
	}

	T const &
	data() const {
		return data_;
	}

	const_iterator
	begin() const {
		return data_.begin();
	}

	const_iterator
	end() const {
		return data_.end();
	}

	iterator
	begin() {
		return data_.begin();
	}

	iterator
	end() {
		return data_.end();
	}

	Size
	size() const {
		return data_.size();
	}

	/// @brief Strict ordering with preference given to the values closest
	/// to the containers begin() element.
	bool
	operator < ( OrderedTuple< T > const & rhs ) const {
		debug_assert( size() == rhs.size() );

		const_iterator lhs_iter( data_.begin()), lhs_iter_end( data_.end() );
		const_iterator rhs_iter( rhs.data_.begin() );

#ifndef NDEBUG
		const_iterator rhs_iter_end( rhs.data_.end() );
#endif

		while ( lhs_iter != lhs_iter_end ) {
			debug_assert( rhs_iter != rhs_iter_end );
			if ( *lhs_iter == *rhs_iter ) {
				++lhs_iter;
				++rhs_iter;
			} else {
				return *lhs_iter < *rhs_iter;
			}
		}
		return false;
	}

	/// @brief Simple comparison operator for the tuple.  Sweeps from begin() to end();
	bool
	operator == ( OrderedTuple< T > const & rhs ) const {
		debug_assert( size() == rhs.size() );

		const_iterator lhs_iter( data_.begin()), lhs_iter_end( data_.end() );
		const_iterator rhs_iter( rhs.data_.begin() );

#ifndef NDEBUG
		const_iterator rhs_iter_end( rhs.data_.end() );
#endif

		while ( lhs_iter != lhs_iter_end ) {
			debug_assert( rhs_iter != rhs_iter_end );
			if ( *lhs_iter == *rhs_iter ) {
				++lhs_iter;
				++rhs_iter;
			} else {
				return false;
			}
		}
		return true;
	}


private:

	container data_;
};

}

#endif
