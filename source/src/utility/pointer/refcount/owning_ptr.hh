// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/owning_ptr.hh
/// @brief  Shared-ownership intrusive reference-counted smart pointer
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Holds a pointer to an object that owning_ptr owns with other smart
///      pointers using the object's reference counting mechanism.
///  @li The object type can be const.
///  @li Can copy construct and assign from pointer, object, or owning_ptr
///      of an assignable type, including from owning_ptr< T > to
///      owning_ptr< T const >.
///  @li Functions owning_ptr_acquire() and owning_ptr_release() must be defined
///      for all types that will be reference counted.


#ifndef INCLUDED_utility_pointer_refcount_owning_ptr_hh
#define INCLUDED_utility_pointer_refcount_owning_ptr_hh


// Unit headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>

// Package headers
#include <utility/pointer/refcount/owning_ptr.functions.hh>

// Project headers
#include <utility/down_cast.hh>

// C++ headers
#include <cassert>
#include <iosfwd>

#include <boost/type_traits/remove_const.hpp>

namespace utility {
namespace pointer {


/// @brief Shared-ownership intrusive reference-counted smart pointer
template< typename T >
class owning_ptr
{


private: // Friends


	template< typename > friend class owning_ptr; // Friendship across object types


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
	owning_ptr() :
		p_( 0 )
	{}


	/// @brief Copy constructor (implicit)
	inline
	owning_ptr( owning_ptr const & a ) :
		p_( a.p_ )
	{
		if ( p_ ) owning_ptr_acquire( p_ );
	}


	/// @brief Assignable copy constructor (implicit)
	template< typename Assignable >
	inline
	owning_ptr( owning_ptr< Assignable > const & a ) :
		p_( a.p_ )
	{
		if ( p_ ) owning_ptr_acquire( p_ );
	}


	/// @brief Object pointer constructor (implicit)
	inline
	owning_ptr( pointer object_p ) :
		p_( object_p )
	{
		if ( p_ ) owning_ptr_acquire( p_ );
	}


	/// @brief Object constructor
	inline
	explicit
	owning_ptr( reference object ) :
		p_( &object )
	{
		if ( p_ ) owning_ptr_acquire( p_ );
	}


	/// @brief Destructor
	inline
	~owning_ptr()
	{
		if ( p_ ) owning_ptr_release( p_ );
	}


public: // Methods: assignment


	/// @brief Copy assignment
	inline
	owning_ptr &
	operator =( owning_ptr const & a )
	{
		// Safe and faster (on average) to skip self-assignment check
		// Faster than the swap method: owning_ptr( a ).swap( *this )
		pointer const po( p_ );
		p_ = a.p_;
		if ( p_ ) owning_ptr_acquire( p_ ); // Add before release in case p_ == a.p_
		if ( po ) owning_ptr_release( po );
		return *this;
	}


	/// @brief Assignable copy assignment
	template< typename Assignable >
	inline
	owning_ptr &
	operator =( owning_ptr< Assignable > const & a )
	{
		// Faster than the swap method: owning_ptr( a ).swap( *this )
		if ( p_ ) owning_ptr_release( p_ ); // Safe if p_ == a.p_ because ref count >= 2
		p_ = a.p_;
		if ( p_ ) owning_ptr_acquire( p_ );
		return *this;
	}


	/// @brief Object pointer assignment
	inline
	owning_ptr &
	operator =( pointer object_p )
	{
		// Faster than the swap method: owning_ptr( object_p ).swap( *this )
		pointer const po( p_ );
		p_ = object_p;
		if ( p_ ) owning_ptr_acquire( p_ ); // Add before release in case p_ == object_p
		if ( po ) owning_ptr_release( po );
		return *this;
	}


	/// @brief Object assignment
	inline
	owning_ptr &
	operator =( reference object )
	{
		// Faster than the swap method: owning_ptr( &object ).swap( *this )
		pointer const po( p_ );
		p_ = &object;
		owning_ptr_acquire( p_ ); // Add before release in case p_ == &object
		if ( po ) owning_ptr_release( po );
		return *this;
	}


public: // Methods: conversions


	/// @brief bool conversion: points to something?
	/// @note Enables unwanted conversions but work-around has compiler-dependencies
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


	/// @brief Raw pointer
	inline
	pointer
	get() const
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


	/// @brief Release ownership and reset
	inline
	void
	reset_to_null()
	{
		if ( p_ ) {
			pointer const po( p_ );
			p_ = 0;
			owning_ptr_release( po );
		}
	}

	/// @brief Release ownership but leave pointed memory alone
	inline
	void
	relinquish_ownership()
	{
		p_ = 0;
	}

	/// @brief Swap
	inline
	void
	swap( owning_ptr & a )
	{
		pointer const po( p_ );
		p_ = a.p_;
		a.p_ = po;
	}

public: // temp helpers for the transition to std::shared_ptr
	inline
	owning_ptr( access_ptr<T> const & a ) :
		p_( a() )
	{
		if ( p_ ) owning_ptr_acquire( p_ );
	}

	/// @brief Release ownership and reset
	inline
	void
	reset()
	{
		if ( p_ ) {
			pointer const po( p_ );
			p_ = 0;
			owning_ptr_release( po );
		}
	}


private: // Fields
#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const {
			bool a = true;
			if( p_ != 0 ) {
				ar << a;
				ar << *p_;
			} else {
				a = false;
				ar << a;
			}

	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version) {
			bool test;
			ar >> test;
			if( test ) {
				pointer a = new T;
				ar >> *( const_cast< typename boost::remove_const<T>::type * > (a) );
				p_ = a;
				owning_ptr_acquire( a );
			}
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif


	/// @brief Pointer to object
	pointer p_;


}; // owning_ptr


// Free Functions


/// @brief owning_ptr == owning_ptr
template< typename T, typename U >
inline
bool
operator ==( owning_ptr< T > const & a, owning_ptr< U > const & b )
{
	return ( a() == b() );
}


/// @brief owning_ptr == pointer
template< typename T >
inline
bool
operator ==( owning_ptr< T > const & a, T const * const b )
{
	return ( a() == b );
}


/// @brief pointer == owning_ptr
template< typename T >
inline
bool
operator ==( T const * const a, owning_ptr< T > const & b )
{
	return ( a == b() );
}


/// @brief owning_ptr != owning_ptr
template< typename T, typename U >
inline
bool
operator !=( owning_ptr< T > const & a, owning_ptr< U > const & b )
{
	return ( a() != b() );
}


/// @brief owning_ptr != pointer
template< typename T >
inline
bool
operator !=( owning_ptr< T > const & a, T const * const b )
{
	return ( a() != b );
}


/// @brief pointer != owning_ptr
template< typename T >
inline
bool
operator !=( T const * const a, owning_ptr< T > const & b )
{
	return ( a != b() );
}


/// @brief owning_ptr < owning_ptr
template< typename T, typename U >
inline
bool
operator <( owning_ptr< T > const & a, owning_ptr< U > const & b )
{
	return ( a() < b() );
}


/// @brief owning_ptr < pointer
template< typename T >
inline
bool
operator <( owning_ptr< T > const & a, T const * const b )
{
	return ( a() < b );
}


/// @brief pointer < owning_ptr
template< typename T >
inline
bool
operator <( T const * const a, owning_ptr< T > const & b )
{
	return ( a < b() );
}


/// @brief Swap two owning_ptr objects
template< typename T >
inline
void
swap( owning_ptr< T > & a, owning_ptr< T > & b )
{
	a.swap( b );
}


/// @brief Get pointer of owning_ptr: needed by boost::mem_fn
template< typename T >
inline
T *
get_pointer( owning_ptr< T > const & a )
{
	return a();
}


/// @brief Static cast an owning_ptr
template< typename T, typename U >
inline
owning_ptr< T >
static_pointer_cast( owning_ptr< U > const & a )
{
	return static_cast< T * >( a() );
}


/// @brief Const cast an owning_ptr
template< typename T, typename U >
inline
owning_ptr< T >
const_pointer_cast( owning_ptr< U > const & a )
{
	return const_cast< T * >( a() );
}


/// @brief Dynamic cast an owning_ptr
template< typename T, typename U >
inline
owning_ptr< T >
dynamic_pointer_cast( owning_ptr< U > const & a )
{
	return dynamic_cast< T * >( a() );
}


/// @brief Down cast an owning_ptr
template< typename T, typename U >
inline
owning_ptr< T >
down_pointer_cast( owning_ptr< U > const & a )
{
	return utility::down_cast< T * >( a() );
}


/// @brief Stream output
template< typename CharT, typename CharTraits, typename T >
inline
std::basic_ostream< CharT, CharTraits > &
operator <<( std::basic_ostream< CharT, CharTraits > & os, owning_ptr< T > const & a )
{
	os << a();
	return os;
}


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_refcount_owning_ptr_HH
