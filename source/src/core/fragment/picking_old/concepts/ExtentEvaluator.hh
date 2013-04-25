// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/picking_old/concepts/ExtentEvaluator.fwd.hh
/// @brief  class demonstrating ExtentEvaluator concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_ExtentEvaluator_hh
#define INCLUDED_core_fragment_picking_old_concepts_ExtentEvaluator_hh


// unit headers
#include <core/fragment/picking_old/concepts/ExtentEvaluator.fwd.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief class demonstrating ExtentEvaluator concept
/// @remarks For demonstration only, do not derive from this class!
template< typename Bookmark, typename Ext >
class ExtentEvaluator {


public:


	typedef Ext Extent;


	/// @brief evaluate an extent of pages
	/// @param[in] extent_begin iterator pointing to the beginning of the Page extent
	/// @param[in] extent_end iterator pointing past the end of the Page extent
	/// @param[out] mark Bookmark that will have evaluation information filled in
	/// @return true if extent is allowed and has been scored, otherwise false
	/// @remarks implement this to get custom extent evaluation and exclusion behavior
	virtual
	bool operator ()( Extent const & extent, Bookmark & mark );


	/// @brief clone this object
	virtual
	ExtentEvaluator * clone();


};


} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragments_picking_old_concepts_ExtentEvaluator_HH */
