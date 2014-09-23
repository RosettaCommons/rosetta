// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCountMI.hh
/// @brief  Interface class for reference-counted multiple inheritance polymorphic classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Intended for use as a base class of polymorphic classes:
///      A template-based approach without virtual destructor is preferred
///      for non-polymorphic classes for efficiency.
///  @li Count value made mutable and reference add/subtract functions made const
///      so that const objects can be held by shared ownership smart pointers.
///  @li Use instead of ReferenceCount if multiple inheritance is being used and a pure
///      interface is desired at the top of the hierarchy.


#ifndef INCLUDED_utility_pointer_ReferenceCountMI_hh
#define INCLUDED_utility_pointer_ReferenceCountMI_hh

#include "platform/types.hh"

// Unit headers
#include <utility/pointer/refcount/ReferenceCountMI.fwd.hh>

// C++ headers
#include <cstddef>


namespace utility {
namespace pointer {


/// @brief Interface class for reference counted polymorphic classes
class ReferenceCountMI
{


private: // Friends


	template< typename T > friend void owning_ptr_acquire( T * );
	template< typename T > friend void owning_ptr_release( T * );


public: // Types


	// Project style
	typedef platform::Size Size;

	// STL/boost style
	typedef  platform::Size size_type;


protected: // Creation


	/// @brief Default constructor
	inline
	ReferenceCountMI()
	{}


	/// @brief Copy constructor
	inline
	ReferenceCountMI( ReferenceCountMI const & )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~ReferenceCountMI()
	{}


protected: // Assignment


	/// @brief Copy assignment
	inline
	ReferenceCountMI &
	operator =( ReferenceCountMI const & )
	{
		return *this;
	}


protected: // Methods


	/// @brief Add a reference: Increment the count
	virtual
	void
	add_ref() const = 0;


	/// @brief Remove a reference: Decrement the count: Self-destruct if count hits zero
	virtual
	void
	remove_ref() const = 0;


public: // Properties


	/// @brief Reference count
	virtual
	Size
	ref_count() const = 0;


}; // ReferenceCountMI


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_ReferenceCountMI_HH
