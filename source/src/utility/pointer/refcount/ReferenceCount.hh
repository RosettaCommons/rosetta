// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCount.hh
/// @brief  Base class for reference-counted single inheritance polymorphic classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Intended for use as a base class of polymorphic classes:
///      A template-based approach without virtual destructor is preferred
///      for non-polymorphic classes for efficiency.
///  @li Count value made mutable and reference add/subtract functions made const
///      so that const objects can be held by shared ownership smart pointers.
///  @li Prefer to ReferenceCountMI if single inheritance is being used and a pure
///      interface is not needed at the top of the hierarchy since the counter
///      update functions are non-virtual here.


#ifndef INCLUDED_utility_pointer_refcount_ReferenceCount_hh
#define INCLUDED_utility_pointer_refcount_ReferenceCount_hh

#ifdef PTR_REFCOUNT

#include <platform/types.hh>

// Unit headers
#include <utility/pointer/ReferenceCount.fwd.hh>

// C++ headers
#include <cassert>
#include <cstddef>

#ifdef MULTI_THREADED
// Boost headers
#include <boost/detail/atomic_count.hpp>
#endif

namespace utility {
namespace pointer {


/// @brief Base class for reference-counted polymorphic classes
class ReferenceCount
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
	ReferenceCount() :
		count_( 0 )
	{}


	/// @brief Copy constructor
	inline
	ReferenceCount( ReferenceCount const & ) :
		count_( 0 ) // New object has no references
	{}


public: // Creation

	void ctor();

	/// @brief Destructor
	inline
	virtual
	~ReferenceCount()
	{
		assert( count_ == 0 ); // Check for dangling references
	}


protected: // Assignment


	/// @brief Copy assignment
	inline
	ReferenceCount &
	operator =( ReferenceCount const & )
	{
		// Assignment doesn't change reference count
		return *this;
	}


private: // Methods


	/// @brief Add a reference: Increment the count
	inline
	void
	add_ref() const
	{
		assert( count_ < max_count_ );
		++count_;
	}


	/// @brief Remove a reference: Decrement the count: Self-destruct if count hits zero
	inline
	void
	remove_ref() const
	{
		assert( count_ > 0 );
		--count_;
		if ( count_ == 0 ) delete this;
	}


public: // Properties


	/// @brief Reference count
	inline
	Size
	ref_count() const
	{
		return count_;
	}


private: // Fields

#ifdef MULTI_THREADED
	/// @brief Max count
	static long const max_count_; // boost uses long.

	/// @brief Reference count
	/// @note Identity semantics: Not copied or assigned
	mutable boost::detail::atomic_count count_;
#else

	/// @brief Max count
	static Size const max_count_;

	/// @brief Reference count
	/// @note Identity semantics: Not copied or assigned
	mutable Size count_;
#endif

}; // ReferenceCount


} // namespace pointer
} // namespace utility


#endif // PTR_REFCOUNT

#endif // INCLUDED_utility_pointer_refcount_ReferenceCount_HH
