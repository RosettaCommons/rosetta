// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/in_place_list.hh
/// @brief  a doubly linked list where elements are shuffled in place.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_in_place_list_hh
#define INCLUDED_utility_in_place_list_hh

// Unit headers
#include <utility/in_place_list.fwd.hh>

// Platform headers
#include <platform/types.hh>

// ObjexxFCL Headers
#include <utility/vector1.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace utility {

template < class T >
class list_element
{
public:
	typedef platform::Size Size;

public:
	list_element() :
		in_list_( false ),
		prev_( 0 ),
		next_( 0 )
	{}

	list_element( T const & val ) :
		in_list_( false ),
		prev_( 0 ),
		next_( 0 ),
		data_( val )
	{}

	inline T const & data() const { return data_; }
	inline T & data() { return data_; }

	inline bool in_list() const { return in_list_; }
	inline Size next() const { return next_; }
	inline Size prev() const { return prev_; }

private:
	void clear() {
		in_list_ = false;
		prev_ = 0;
		next_ = 0;
	}

public:
	friend class in_place_list< T >;

private:
	bool in_list_;
	Size prev_;
	Size next_;
	T    data_;
};

template < class T >
class in_place_list : public utility::pointer::ReferenceCount {
public:
	typedef platform::Size Size;

public:
	// @brief default constructor; requires a resize
	in_place_list() : head_( 0 ), tail_( 0 ) {}

	// @brief Size constructor; list can be used straight away, unlike with the default constructor
	in_place_list( Size n_items ) : elements_( n_items ), head_( 0 ), tail_( 0 ) {}

	/// @brief Size + value constructor; list can be used straight away, unlike with the default constructor
	in_place_list( Size n_items, T const & val ) : elements_( n_items, val ), head_( 0 ), tail_( 0 ) {}

	Size
	size() const {
		return elements_.size();
	}

	/// @brief Clear the list entirely, and resize to be larger
	void
	resize( Size n_items ) {
		elements_.resize( n_items );
		for ( Size ii = 1; ii <= n_items; ++ii ) elements_[ ii ].clear();
		head_ = 0; tail_ = 0;
	}

	/// @brief O(k) clear -- iterates from head to end, clears each element in the list
	void
	clear() {
		for ( Size ii = head_, iiend = end(); ii != iiend; /* no increment */ ) {
			Size iinext = elements_[ ii ].next();
			elements_[ ii ].clear();
			ii = iinext;
		}
		head_ = 0;
		tail_ = 0;
	}

	inline
	list_element< T > const &
	operator[] ( Size index ) const {
		return elements_[ index ];
	}

	inline
	list_element< T > &
	operator[] ( Size index ) {
		return elements_[ index ];
	}

	inline
	void move_to_front( Size index ) {
		if ( head_ == index ) return;
		if ( elements_[ index ].in_list() ) {
			extract( index );
		}
		set_head( index );
	}

	inline
	void move_to_back( Size index ) {
		if ( tail_ == index ) return;
		if ( elements_[ index ].in_list() ) {
			extract( index );
		}
		set_tail( index );

	}

	inline
	void remove( Size index ) {
	debug_assert( elements_[ index ].in_list() );
		extract( index );
		elements_[ index ].prev_ = 0;
		elements_[ index ].next_ = 0;
		elements_[ index ].in_list_ = false;
	}

	inline Size head() const { return head_; }
	inline Size tail() const { return tail_; }
	inline Size end() const { return 0; } // iterate from head to end, or from tail to end; either way, end at end()

private:

	/// @brief Remove an element from its place in the list; update the prev_ and next_
	/// pointers of its surrounding elements.  Also update the head_ and tail_
	/// pointers for the list, if necessary
	void extract( Size index ) {
	debug_assert( elements_[ index ].in_list() );
		Size prev = elements_[ index ].prev();
		Size next = elements_[ index ].next();
		if ( prev != 0 ) {
		debug_assert( elements_[ prev ].next_ == index );
			elements_[ prev ].next_ = next;
		} else {
		debug_assert( head_ == index );
			head_ = next;
		}
		if ( next != 0 ) {
		debug_assert( elements_[ next ].prev_ == index );
			elements_[ next ].prev_ = prev;
		} else {
		debug_assert( tail_ == index );
			tail_ = prev;
		}
	}

	/// @brief Simply set the head-element in the list.  Does not clean up
	/// pointers for the surrounding elements if the requested element is
	/// already a member of the list.
	void set_head( Size index ) {
		if ( head_ != 0 ) {
		debug_assert( elements_[ head_ ].prev_ == 0 );
			elements_[ head_ ].prev_ = index;
		}
		elements_[ index ].next_ = head_;
		elements_[ index ].prev_ = 0;
		elements_[ index ].in_list_ = true;
		head_ = index;
		if ( tail_ == 0 ) {
			tail_ = head_;
		}

	}

	/// @brief Simply set the tail-element in the list.  Does not clean up
	/// pointers for the surrounding elements if the requested element is
	/// already a member of the list.
	void set_tail( Size index ) {
		if ( tail_ != 0 ) {
		debug_assert( elements_[ tail_ ].next_ == 0 );
			elements_[ tail_ ].next_ = index;
		}
		elements_[ index ].prev_ = tail_;
		elements_[ index ].next_ = 0;
		elements_[ index ].in_list_ = true;
		tail_ = index;
		if ( head_ == 0 ) {
			head_ = tail_;
		}
	}

private:
	typename utility::vector1< list_element< T > > elements_;
	Size head_;
	Size tail_;

}; // class heap

}

#endif
