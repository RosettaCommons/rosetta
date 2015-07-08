// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/frags/heap.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1A.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end

// #include <ObjexxFCL/FArray2D.hh>
// #include <ObjexxFCL/FArray3D.hh>
// #include <ObjexxFCL/FArray4D.hh>
//#include <ObjexxFCL/format.hh>

// C++ Headers
//#include <algorithm>
//#include <iostream>

namespace protocols {
namespace frags {

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//  heap functions
//     heap is a list integers the value associated with each integer is
//     in the coheap (or covalue) array. Comparisons <=> are made between
//     the covalues.  The heap entries are thus sorted on the basis of
//     their covalues.  In typical usage, this is an index type sort
//     the heap values are the indicies, and the covalues the value.
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
///
/// @brief
/// sets up an empty heap and stores the dimensioned size
///
/// @details
///
/// @param  heap - [in/out]? -
/// @param  coheap - [in/out]? -
/// @param  max_items - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_init(
	FArray1A_int heap,
	FArray1A_float coheap,
	int max_items
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );

	heap(-1) = 0;
	heap(-2) = max_items;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// modifes heap and last_val
/// return val and err.
///
/// @details
///
/// @param  heap - [in/out]? - convert to zero offset matrix
/// @param  coheap - [in/out]? -
/// @param  val - [in/out]? -
/// @param  coval - [in/out]? -
/// @param  err - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_extract(
	FArray1A_int heap, // convert to zero offset matrix
	FArray1A_float coheap,
	int & val,
	float & coval,
	bool & err
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );

	err = true;
	if ( heap(-1) < 1 ) return;

	int temp_val = heap(0);
	float temp_coval = coheap(0);
	err = false;
	--heap(-1);

	if ( heap(-1) == 0 ) { // special case for single value in heap
		val = temp_val; // PB bugfix
		coval = temp_coval;
		return;
	}

	heap(0) = heap(heap(-1)); // move last value to front
	coheap(0) = coheap(heap(-1));

	heap_down(heap,coheap,1);
//  we use a temporary copy so that this can be done in place if need be
	val = temp_val;
	coval = temp_coval;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// modifes heap and last_dummy, inserts val, returns err
/// requires heap_max to be previously set via heap_init
///
/// @details
///
/// @param  heap - [in/out]? - convert to zero offset matrix
/// @param  coheap - [in/out]? -
/// @param  val - [in/out]? -
/// @param  coval - [in/out]? -
/// @param  err - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_insert(
	FArray1A_int heap, // convert to zero offset matrix
	FArray1A_float coheap,
	int val,
	float coval,
	bool & err
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );

	if ( heap(-1) >= heap(-2) ) { // list is full, use replace instead
		err = true;
		if ( coheap(0) < coval ) heap_replace(heap,coheap,val,coval);
		return;
	}

	err = false;
	heap(heap(-1)) = val; // empty spot on end (zero offset)
	coheap(heap(-1)) = coval;

	++heap(-1);
	heap_up(heap,coheap,heap(-1));

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  heap - [in/out]? - convert to zero offset matrix
/// @param  coheap - [in/out]? -
/// @param  val - [in/out]? -
/// @param  coval - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_replace(
	FArray1A_int heap, // convert to zero offset matrix
	FArray1A_float coheap,
	int val,
	float coval
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );

// modifes heap

	//bool err;
	//err = false;  set but never used ~Labonte

	heap(0) = val; // overwrite the lowest element
	coheap(0) = coval;
	heap_down(heap,coheap,1);
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  heap - [in/out]? - convert to zero offset matrix
/// @param  coheap - [in/out]? -
/// @param  index_in - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_down(
	FArray1A_int heap, // convert to zero offset matrix
	FArray1A_float coheap,
	int index_in
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );
	float coiv,/*cocv,*/cocv2;
	int indx,/*child,*/iv,/*cv,*/cv2,last;
	indx = index_in-1; // convert to zero offset matrix
	last = heap(-1)-1; // convert to zero offset matrix

	if ( last <= 0 ) return; // empty or single element
	if ( indx > last ) return; // dumbass

	iv = heap(indx); // the inserted value
	coiv = coheap(indx);

	while ( indx < last ) {
		int child = 2*indx+1;

		if ( child > last ) break;

		int cv  = heap(child);
		float cocv = coheap(child);

		if ( child < last ) {
			cv2 = heap (child+1);
			cocv2 = coheap(child+1);

			if ( cocv2 < cocv ) {
				cv = cv2;
				cocv = cocv2;

				++child;
			}
		}

		if ( coiv <= cocv ) break;
		coheap(indx) = cocv;
		heap(indx) = cv;
		indx = child;
	}

	heap(indx) = iv;
	coheap(indx) = coiv;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  heap - [in/out]? - convert to zero offset matrix
/// @param  coheap - [in/out]? -
/// @param  index_in - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
heap_up(
	FArray1A_int heap, // convert to zero offset matrix
	FArray1A_float coheap,
	int & index_in
)
{
	heap.dimension( SRange( -2, star ) );
	coheap.dimension( SRange( 0, star ) );
	float covalue;//,copv;

	int indx,/*parent,*/value;//,pv;
	indx = index_in-1; // convert to zero offset matrix


	value = heap(indx);
	covalue = coheap(indx);

	while ( indx != 0 ) {
		int parent = static_cast< int >((indx-1)/2);
		int pv = heap(parent);
		float copv = coheap(parent);
		if ( copv < covalue ) break;
		coheap(indx) = copv;
		heap(indx) = pv;
		indx = parent;
	}

	coheap(indx) = covalue;
	heap(indx) = value;
}


} // ns frags
} // ns protocols
