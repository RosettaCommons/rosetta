// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file heap.cc
/// @brief class definition for a heap object based on Charlie
/// Strauss's heap code ported over from rosetta++.
/// @author James Thompson
/// @author Charlie Strauss

// Rosetta Headers
#include <utility/heap.hh>

// ObjexxFCL Headers

// C++ Headers
namespace utility {

// @brief Auto-generated virtual destructor
heap::~heap() {}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//  heap functions
//     heap is a list integers the value associated with each integer is
//     in the coheap (or covalue) array. Comparisons <=> are made between
//     the covalues.  The heap entries are thus sorted on the basis of
//     their covalues.  In typical usage, this is an index type sort
//     the heap values are the indices, and the covalues the value.
//
//     heap is simply an array with two extra storage elments at the
//     front the first is the max dimension of the heap, the second is
//     the current number of entries (the third is the minimum value in
//     heap and the start of the heap). When dimensioning space for it be
//     sure to add 2 elements more than you think you need.
//
//     heap_init    set up an empty heap
//     heap_insert  insert a new value into the heap
//     heap_extract extract the lowset value (always located in heap(3))
//
//     heap_replace replace the lowest value
//        (equivalent to heap_extract; heap_insert  but faster)
//        If you call heap_insert with a full heap (ie last = maxsize) then
//        heap_replace gets called instead.

//     charlie strauss 1999
//------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// @brief sets up an empty heap and stores the dimensioned size
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_init(
	int max_items
)
{
	heap_size() = 0;
	heap_capacity() = max_items;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// modifes heap and last_val return val and err.
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_extract(
	int & val,
	float & coval,
	bool & err
)
{
	err = true;
	//if ( heap_(-1) < 1 ) return;
	if ( heap_size() < 1 ) return;

	//int temp_val = heap_(0);
	//float temp_coval = coheap_(0);
	int temp_val = heap_[0];
	float temp_coval = coheap_[0];
	err = false;
	--heap_size();

	if ( heap_size() == 0 ) { // special case for single value in heap
		val = temp_val; // PB bugfix
		coval = temp_coval;
		return;
	}

	//heap_(0) = heap_(heap_(-1)); // move last value to front
	//coheap_(0) = coheap_(heap_(-1));
	heap_  [0] = heap_  [ heap_size() ]; // move last value to front
	coheap_[0] = coheap_[ heap_size() ];

	heap_down(1);
	//  we use a temporary copy so that this can be done in place if need be
	val = temp_val;
	coval = temp_coval;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// modifies heap and last_dummy, inserts val, returns err
/// requires heap_max to be previously set via heap_init
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_insert(
	int val,
	float coval,
	bool & err
)
{
	if ( heap_size() >= heap_capacity() ) {
		// list is full, use replace instead
		err = true;
		if ( coheap_[0] < coval ) heap_replace(val,coval);
		return;
	}

	err = false;
	//heap_(heap_(-1)) = val; // empty spot on end (zero offset)
	//coheap_(heap_(-1)) = coval;
	heap_  [ heap_size() ] = val; // empty spot on end (zero offset)
	coheap_[ heap_size() ] = coval;

	++heap_size();
	heap_up( heap_size() );

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_replace(
	int val,
	float coval
)
{
	// modifes heap

	//bool err;
	//err = false;

	heap_[0] = val; // overwrite the lowest element
	coheap_[0] = coval;
	heap_down(1);
}

void
heap::reset_coval( int val, float coval )
{
	int index = index_for_val( val );
	if ( index == -1 ) return; // ERROR

	if ( coheap_[ index ] < coval ) {
		increase_coval( index, coval );
	} else if ( coheap_[ index ] > coval ) {
		decrease_coval( index, coval );
	} // else, noop -- if the old coval == new coval, don't shuffle any part of the heap
}


/// @brief returns the smallest covalue stored in the heap.
float
heap::heap_head() const
{
	return coheap_[ 0 ];
}

float
heap::coval( int index ) const
{
	return coheap_[ index ];
}

int
heap::val( int index ) const
{
	return heap_[ index ];
}

int
heap::size() const
{
	return heap_[ heap_.size() - 1 ];
}

int
heap::capacity() const
{
	return heap_[ heap_.size() - 2 ];
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_down(
	int index_in
)
{
	float coiv,cocv2;
	int indx,iv,/*cv,*/cv2,last;
	indx = index_in-1; // convert to zero offset matrix
	last = heap_size() - 1; // convert to zero offset matrix

	if ( last <= 0 ) return; // empty or single element
	if ( indx > last ) return; // dumbass

	iv   = heap_  [indx]; // the inserted value
	coiv = coheap_[indx];

	while ( indx < last ) {
		int child = 2*indx+1;

		if ( child > last ) break; // GHAA goto L20; // loop escape

		int cv  = heap_[child];
		float cocv = coheap_[child];

		if ( child < last ) {
			cv2 = heap_[child+1];
			cocv2 = coheap_[child+1];

			if ( cocv2 < cocv ) {
				cv = cv2;
				cocv = cocv2;

				++child;
			}
		}

		if ( coiv <= cocv ) break; /// GHAA goto L20; // loop escape
		coheap_[indx] = cocv;
		heap_[indx] = cv;
		indx = child;
	}

	/// WTF? L20:; // loop escape
	heap_[indx] = iv;
	coheap_[indx] = coiv;
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void
heap::heap_up(
	int index_in
)
{
	float covalue;//,copv;

	int indx,value;
	indx = index_in-1; // convert to zero offset matrix


	value = heap_[indx];
	covalue = coheap_[indx];

	while ( indx != 0 ) {
		int parent = static_cast< int >((indx-1)/2);
		int pv = heap_[parent];
		float copv = coheap_[parent];
		if ( copv < covalue ) break; // GHAA goto L20; // loop escape
		coheap_[indx] = copv;
		heap_[indx] = pv;
		indx = parent;
	}

	// GHAA L20:; // loop escape
	coheap_[indx] = covalue;
	heap_[indx] = value;
}

int &
heap::heap_size()
{
	return heap_[ heap_.size() - 1 ];
}

int &
heap::heap_capacity()
{
	return heap_[ heap_.size() - 2 ];
}


void
heap::decrease_coval( int index, float coval )
{
	coheap_[ index ] = coval;
	heap_up( ++index ); // so the index can be decremented again...
}

void
heap::increase_coval( int index, float coval )
{
	coheap_[ index ] = coval;
	heap_down( ++index ); // so the index can be decremented again...
}

int
heap::index_for_val( int val ) {
	for ( int ii = 0, ii_end = heap_size(); ii < ii_end; ++ii ) {
		if ( heap_[ ii ] == val ) {
			return ii;
		}
	}
	return -1;
}


} // ns utility
