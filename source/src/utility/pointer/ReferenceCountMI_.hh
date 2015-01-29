// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCountMI_.hh
/// @brief  Base class for reference-counted multiple inheritance polymorphic classes
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


#ifndef INCLUDED_utility_pointer_ReferenceCountMI__hh
#define INCLUDED_utility_pointer_ReferenceCountMI__hh


// Unit headers
#include <utility/pointer/ReferenceCountMI.hh>

// C++ headers
#include <utility/assert.hh>


namespace utility {
namespace pointer {


/// @brief Base class for reference counted polymorphic classes
class ReferenceCountMI_ :
	virtual public ReferenceCountMI
{


private: // Friends


	template< typename T > friend void owning_ptr_acquire( T * );
	template< typename T > friend void owning_ptr_release( T * );


private: // Types


	typedef  ReferenceCountMI  Super;


protected: // Creation


	/// @brief Default constructor
	inline
	ReferenceCountMI_() :
		count_( 0 )
	{}


	/// @brief Copy constructor
	inline
	ReferenceCountMI_( ReferenceCountMI_ const & ) :
		Super(),
		count_( 0 ) // New object has no references
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~ReferenceCountMI_()
	{
	debug_assert( count_ == 0 ); // Check for dangling references
	}


protected: // Assignment


	/// @brief Copy assignment
	inline
	ReferenceCountMI_ &
	operator =( ReferenceCountMI_ const & )
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
	debug_assert( count_ < max_count_ );
		++count_;
	}


	/// @brief Remove a reference: Decrement the count: Self-destruct if count hits zero
	inline
	void
	remove_ref() const
	{
	debug_assert( count_ > 0 );
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


	/// @brief Max count
	static Size const max_count_;

	/// @brief Reference count
	/// @note Identity semantics: Not copied or assigned
	mutable Size count_;


}; // ReferenceCountMI_


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_ReferenceCountMI__HH
