// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/concepts/Bookmark.hh
/// @brief  struct demonstrating Bookmark concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_Bookmark_hh
#define INCLUDED_core_fragment_picking_old_concepts_Bookmark_hh

// unit headers
#include <core/fragment/picking_old/concepts/Bookmark.fwd.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief struct demonstrating Bookmark concept
/// @remarks For demonstration only, do not derive from this struct.
template< typename PageConstIterator >
struct Bookmark {


public: // comparator


	/// @brief '<' comparison
	bool operator <( Bookmark const & rval ) const;


public: // data


	/// @brief points to the beginning of the fragment
	PageConstIterator extent_begin;


	/// @brief points just beyond the end of the fragment
	PageConstIterator extent_end;


};


} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_concepts_Bookmark_HH */
