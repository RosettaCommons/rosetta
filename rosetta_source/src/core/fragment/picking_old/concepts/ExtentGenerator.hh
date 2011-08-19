// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/concepts/ExtentGenerator.hh
/// @brief class demonstrating ExtentGenerator concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_ExtentGenerator_hh
#define INCLUDED_core_fragment_picking_old_concepts_ExtentGenerator_hh

// unit headers
#include <core/fragment/picking_old/concepts/ExtentGenerator.fwd.hh>

// package headers
#include <core/fragment/picking_old/concepts/Extent.hh>
//#include <core/fragment/picking_old/PageIterator.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {

/*
/// @brief class demonstrating ExtentGenerator concept
/// @remarks for demonstration only, do not derive from this class!
template< typename Ext >
class ExtentGenerator {


public:


	typedef Ext Extent;


	/// @brief given the beginning of an extent, return the desired extent end
	/// @param extent_begin the beginning Page of the extent
	/// @param book_end points just beyond the last Page of the extent
	/// @remarks implement this to get custom extent generation
	Extent operator ()( PageIterator extent_begin, PageIterator book_end );


	/// @brief clone this object
	virtual
	ExtentGenerator * clone();


};
*/

} // concepts
} // picking_old
} // fragment
} // core

#endif /* INCLUDED_core_fragment_picking_old_concepts_ExtentGenerator_HH */
