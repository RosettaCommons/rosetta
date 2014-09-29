// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/std/ReferenceCount.hh
/// @brief  Dummy base class
/// @author


#ifndef INCLUDED_utility_pointer_std_ReferenceCount_hh
#define INCLUDED_utility_pointer_std_ReferenceCount_hh

#ifdef PTR_STD

#include <platform/types.hh>

// Unit headers
#include <utility/pointer/ReferenceCount.fwd.hh>

namespace utility {
namespace pointer {


/// @brief Base class for reference-counted polymorphic classes
class ReferenceCount
{

public: // Types

	// Project style
	typedef platform::Size Size;

	// STL/boost style
	typedef platform::Size size_type;

public: // Creation


	/// @brief Default constructor
	inline
	ReferenceCount()
	{}


	inline
	virtual
	~ReferenceCount()
	{}


}; // ReferenceCount


} // namespace pointer
} // namespace utility

#endif // PTR_STD

#endif // INCLUDED_utility_pointer_std_ReferenceCount_HH
