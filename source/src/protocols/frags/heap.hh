// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


#ifndef INCLUDED_protocols_frags_heap_hh
#define INCLUDED_protocols_frags_heap_hh


// ObjexxFCL Headers

#include <ObjexxFCL/FArray1A.fwd.hh>


// C++ Headers
//#include <iosfwd>

namespace protocols {
namespace frags {

void
heap_init(
	ObjexxFCL::FArray1A_int heap,
	ObjexxFCL::FArray1A_float coheap,
	int max_items
);


void
heap_extract(
	ObjexxFCL::FArray1A_int heap, // convert to zero offset matrix
	ObjexxFCL::FArray1A_float coheap,
	int & val,
	float & coval,
	bool & err
);


void
heap_insert(
	ObjexxFCL::FArray1A_int heap, // convert to zero offset matrix
	ObjexxFCL::FArray1A_float coheap,
	int val,
	float coval,
	bool & err
);


void
heap_replace(
	ObjexxFCL::FArray1A_int heap, // convert to zero offset matrix
	ObjexxFCL::FArray1A_float coheap,
	int val,
	float coval
);


void
heap_down(
	ObjexxFCL::FArray1A_int heap, // convert to zero offset matrix
	ObjexxFCL::FArray1A_float coheap,
	int index_in
);


void
heap_up(
	ObjexxFCL::FArray1A_int heap, // convert to zero offset matrix
	ObjexxFCL::FArray1A_float coheap,
	int & index_in
);


} // ns frags
} // ns protocols


#endif
