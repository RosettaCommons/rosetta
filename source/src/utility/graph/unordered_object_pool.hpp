// Copyright (C) 2000, 2001 Stephen Cleary
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org for updates, documentation, and revision history.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef BOOST_UNORDERED_OBJECT_POOL_HPP
#define BOOST_UNORDERED_OBJECT_POOL_HPP

#include <utility/graph/unordered_object_pool.fwd.hpp>

#ifdef PYROSETTA
#include <boost/pool/pool.hpp>
#endif

// boost::pool

// The following code will be put into Boost.Config in a later revision
#if defined(BOOST_MSVC) || defined(__KCC)
# define BOOST_NO_TEMPLATE_CV_REF_OVERLOADS
#endif

// The following code might be put into some Boost.Config header in a later revision
#ifdef __BORLANDC__
# pragma option push -w-inl
#endif

// There are a few places in this file where the expression "this->m" is used.
// This expression is used to force instantiation-time name lookup, which I am
//   informed is required for strict Standard compliance.  It's only necessary
//   if "m" is a member of a base class that is dependent on a template
//   parameter.
// Thanks to Jens Maurer for pointing this out!

// APL This class differs from the object-pool class in that its destructor does
// not invoke the destructors of un-freed objects that are still living in the pool.
// For destructors to be called, objects must be explicitly destroyed with a call
// to destroy().  The trade-off is hassle for speed.  In code where lots of small objects
// are being created and destroyed on a regular basis, this class is dramatically
// faster.

namespace boost {

// T must have a non-throwing destructor
template <typename T, typename UserAllocator>
class unordered_object_pool: protected pool<UserAllocator>
{
public:
	typedef T element_type;
	typedef UserAllocator user_allocator;
	typedef typename pool<UserAllocator>::size_type size_type;
	typedef typename pool<UserAllocator>::difference_type difference_type;

protected:
	pool<UserAllocator> & store() { return *this; }
	const pool<UserAllocator> & store() const { return *this; }

	// for the sake of code readability :)
	static void * & nextof(void * const ptr)
	{ return *(static_cast<void **>(ptr)); }

public:
	// This constructor parameter is an extension!
	explicit unordered_object_pool(const size_type next_size = 32)
	:pool<UserAllocator>(sizeof(T), next_size) { }

	~unordered_object_pool();

	// Returns 0 if out-of-memory
	element_type * malloc()
	{ return static_cast<element_type *>(store().malloc()); }
	void free(element_type * const chunk)
	{ store().free(chunk); }
	bool is_from(element_type * const chunk) const
	{ return store().is_from(chunk); }

	element_type * construct()
	{
		element_type * const ret = malloc();
		if ( ret == 0 ) {
			return ret;
		}
		/*try*/ { new (ret) element_type(); }
		/*catch (...) { free(ret); throw; }*/
		return ret;
	}

	// Include automatically-generated file for family of template construct()
	//  functions
#ifndef BOOST_NO_TEMPLATE_CV_REF_OVERLOADS
	#   include <boost/pool/detail/pool_construct.ipp>
#else
#   include <boost/pool/detail/pool_construct_simple.ipp>
#endif

	void destroy(element_type * const chunk)
	{
		chunk->~T();
		free(chunk);
	}

	// These functions are extensions!
	size_type get_next_size() const { return store().get_next_size(); }
	void set_next_size(const size_type x) { store().set_next_size(x); }
};

template <typename T, typename UserAllocator>
unordered_object_pool<T, UserAllocator>::~unordered_object_pool()
{
}

} // namespace boost

// The following code might be put into some Boost.Config header in a later revision
#ifdef __BORLANDC__
# pragma option pop
#endif

#endif
