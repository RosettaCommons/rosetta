// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/concepts/Extent.hh
/// @brief  class demonstrating the Extent concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_Extent_hh
#define INCLUDED_core_fragment_picking_old_concepts_Extent_hh

// unit headers
#include <core/fragment/picking_old/concepts/Extent.fwd.hh>

// type headers
#include <core/types.hh>

// C++ headers
#include <iterator>


namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief class demonstrating the Extent concept
/// @remarks Class is usable as a concrete implementation of Extent.
template< typename PageIter >
struct Extent {


public: // typedefs


	typedef core::Size Size;
	typedef PageIter PageIterator;


public: // convenience


	/// @brief compute distance (effectively the length of the extent)
	///  from begin -> end
	inline
	Size distance() const {
		return std::distance( begin, end );
	}


public: // data


	/// @brief points to the beginning of the extent
	PageIterator begin;


	/// @brief points just past the end of the extent
	PageIterator end;


	/// @brief true if extent is to be evaluated by ExtentEvaluators
	///  during Librarian::catalog(), false if it is to be skipped
	/// @remarks ExtentEvaluators do not need to check this flag!
	///  Librarian will only pass valid (true) extents to an evaluator,
	///  invalid (false) extents will be skipped.
	bool valid;


};


} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_concepts_Extent_HH */
