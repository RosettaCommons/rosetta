// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/access_ptr.hh
/// @brief  Non-owning access smart pointer
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Holds a pointer to an object for access but not ownership.
///  @li The object type can be const.
///  @li Can copy construct and assign from pointer, object, or access_ptr
///      of an assignable type, including from access_ptr< T > to
///      access_ptr< T const >.


#ifndef INCLUDED_utility_pointer_refcount_access_ptr_hh
#define INCLUDED_utility_pointer_refcount_access_ptr_hh


// Unit headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>

// Project headers
#include <utility/down_cast.hh>

// C++ headers
#include <cassert>
#include <iosfwd>


namespace utility {
namespace pointer {


/// @brief Non-owning access smart pointer
template< typename T >
class access_ptr
{


private: // Friends


	template< typename > friend class access_ptr; // Friendship across object types


public: // Types


	// Project style
	typedef  T          Value;
	typedef  T &        Reference;
	typedef  T const &  ConstReference;
	typedef  T *        Pointer;
	typedef  T const *  ConstPointer;

	// STL/boost style
	typedef  T          value_type;
	typedef  T          element_type; // For boost::intrusive_ptr compatibility
	typedef  T &        reference;
	typedef  T const &  const_reference;
	typedef  T *        pointer;
	typedef  T const *  const_pointer;


public: // Creation


	/// @brief Default constructor
	inline
	access_ptr() :
		p_( 0 )
	{}


	/// @brief Copy constructor (implicit)
	inline
	access_ptr( access_ptr const & a ) :
		p_( a.p_ )
	{}


	/// @brief Assignable copy constructor (implicit)
	template< typename U >
	inline
	access_ptr( access_ptr< U > const & a ) :
		p_( a.p_ )
	{}


	/// @brief Object pointer constructor (implicit)
	inline
	access_ptr( pointer object_p ) :
		p_( object_p )
	{}


	/// @brief Object constructor
	inline
	explicit
	access_ptr( reference object ) :
		p_( &object )
	{}


	/// @brief Destructor: Non-owning => doesn't delete object
	inline
	~access_ptr()
	{}


public: // Methods: assignment


	/// @brief Copy assignment
	inline
	access_ptr &
	operator =( access_ptr const & a )
	{
		// Safe and faster (on average) to skip self-assignment check
		p_ = a.p_;
		return *this;
	}


	/// @brief Assignable copy assignment
	template< typename U >
	inline
	access_ptr &
	operator =( access_ptr< U > const & a )
	{
		p_ = a.p_;
		return *this;
	}


	/// @brief Object pointer assignment
	inline
	access_ptr &
	operator =( pointer object_p )
	{
		p_ = object_p;
		return *this;
	}


	/// @brief Object assignment
	inline
	access_ptr &
	operator =( reference object )
	{
		p_ = &object;
		return *this;
	}


public: // Methods: conversions


	/// @brief bool conversion: points to something?
	/// @note  Enables unwanted conversions but work-around has compiler-dependencies
	inline
	operator bool() const
	{
		return ( p_ != 0 );
	}


public: // Methods: operators


	/// @brief Dereference
	inline
	reference
	operator *() const
	{
		assert( p_ != 0 );
		return *p_;
	}


	/// @brief Indirection
	inline
	pointer
	operator ->() const
	{
		assert( p_ != 0 );
		return p_;
	}


	/// @brief Raw pointer
	inline
	pointer
	operator ()() const
	{
		return p_;
	}


	/// @brief Points to nothing? (some compilers need this)
	inline
	bool
	operator !() const
	{
		return ( p_ == 0 );
	}


public: // Methods


	/// @brief Raw pointer
	inline
	pointer
	get() const
	{
		return p_;
	}


	/// @brief Reset
	inline
	void
	reset_to_null()
	{
		p_ = 0;
	}


	/// @brief Swap
	inline
	void
	swap( access_ptr & a )
	{
		pointer const po( p_ );
		p_ = a.p_;
		a.p_ = po;
	}


public: // temp helpers for the transition to std::weak_ptr
	inline
	bool expired() const {
		return !p_;
	}

	inline
	void reset() {
		p_ = 0;
	}

	inline
	access_ptr( owning_ptr< T > const & a ) :
		p_( a() )
	{
	}

	inline
	owning_ptr<T> lock() const
	{
		return owning_ptr<T>(*this);
	}

private: // Fields


	/// @brief Pointer to object
	pointer p_;


}; // access_ptr


// Free Functions


/// @brief access_ptr == access_ptr
template< typename T, typename U >
inline
bool
operator ==( access_ptr< T > const & a, access_ptr< U > const & b )
{
	return ( a() == b() );
}


/// @brief access_ptr == pointer
template< typename T >
inline
bool
operator ==( access_ptr< T > const & a, T const * const b )
{
	return ( a() == b );
}


/// @brief pointer == access_ptr
template< typename T >
inline
bool
operator ==( T const * const a, access_ptr< T > const & b )
{
	return ( a == b() );
}


/// @brief access_ptr != access_ptr
template< typename T, typename U >
inline
bool
operator !=( access_ptr< T > const & a, access_ptr< U > const & b )
{
	return ( a() != b() );
}


/// @brief access_ptr != pointer
template< typename T >
inline
bool
operator !=( access_ptr< T > const & a, T const * const b )
{
	return ( a() != b );
}


/// @brief pointer != access_ptr
template< typename T >
inline
bool
operator !=( T const * const a, access_ptr< T > const & b )
{
	return ( a != b() );
}


/// @brief access_ptr < access_ptr
template< typename T, typename U >
inline
bool
operator <( access_ptr< T > const & a, access_ptr< U > const & b )
{
	return ( a() < b() );
}


/// @brief access_ptr < pointer
template< typename T >
inline
bool
operator <( access_ptr< T > const & a, T const * const b )
{
	return ( a() < b );
}


/// @brief pointer < access_ptr
template< typename T >
inline
bool
operator <( T const * const a, access_ptr< T > const & b )
{
	return ( a < b() );
}


/// @brief Swap two access_ptr objects
template< typename T >
inline
void
swap( access_ptr< T > & a, access_ptr< T > & b )
{
	a.swap( b );
}


/// @brief Get pointer of access_ptr: needed by boost::mem_fn
template< typename T >
inline
T *
get_pointer( access_ptr< T > const & a )
{
	return a();
}


/// @brief Static cast an access_ptr
template< typename T, typename U >
inline
access_ptr< T >
static_pointer_cast( access_ptr< U > const & a )
{
	return static_cast< T * >( a() );
}


/// @brief Const cast an access_ptr
template< typename T, typename U >
inline
access_ptr< T >
const_pointer_cast( access_ptr< U > const & a )
{
	return const_cast< T * >( a() );
}


/// @brief Dynamic cast an access_ptr
template< typename T, typename U >
inline
access_ptr< T >
dynamic_pointer_cast( access_ptr< U > const & a )
{
	return dynamic_cast< T * >( a() );
}


/// @brief Down cast an access_ptr
template< typename T, typename U >
inline
access_ptr< T >
down_pointer_cast( access_ptr< U > const & a )
{
	return utility::down_cast< T * >( a() );
}


/// @brief Stream output
template< typename CharT, typename CharTraits, typename T >
inline
std::basic_ostream< CharT, CharTraits > &
operator <<( std::basic_ostream< CharT, CharTraits > & os, access_ptr< T > const & a )
{
	os << a();
	return os;
}

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( access_ptr< T > const & a, access_ptr< U > const & b )
{
	return a && b && a() == b();
}

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( access_ptr< T > & a, owning_ptr< U > const & b )
{
	return a && a() == b();
}

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( access_ptr< T > & a,  U* const b )
{
	return a && a() == b;
}

} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_refcount_access_ptr_HH
