///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) 2010-2011  Artyom Beilis (Tonkikh) <artyomtnk@yahoo.com>
//
//  Distributed under:
//
//                   the Boost Software License, Version 1.0.
//              (See accompanying file LICENSE_1_0.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt)
//
//  or (at your opinion) under:
//
//                               The MIT License
//                 (See accompanying file MIT.txt or a copy at
//              http://www.opensource.org/licenses/mit-license.php)
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CPPDB_REF_PTR_H
#define CPPDB_REF_PTR_H
#include <cppdb/errors.h>
#include <cppdb/atomic_counter.h>

namespace cppdb {
	///
	/// \brief This is a smart intrusive reference counting pointer that throws a error on empty
	/// access.
	///
	/// The T should follow these concepts:
	///
	/// \code
	///   T::add_ref() // increases reference count
	///   IntegerType  T::del_ref() // decreases reference count and returns its current value
	///   static T::dispose(T *) // destroys the object
	/// \endcode
	///
	template<typename T>
	class ref_ptr {
	public:
		///
		/// Default create a new object, if v is not null increases its reference count and stores it
		///
		ref_ptr(T *v=0) : p(0)
		{
			reset(v);
		}
		///
		/// Dereference the object, if reference count goes to 0 destroys it calling T::dispose
		///
		~ref_ptr()
		{
			reset();
		}
		///
		/// Copy a pointer
		///
		ref_ptr(ref_ptr const &other) : p(0)
		{
			reset(other.p);
		}
		///
		/// Assign a pointer
		///
		ref_ptr const &operator=(ref_ptr const &other)
		{
			reset(other.p);
			return *this;
		}
// Borland warns on assignments using operator=(ref_ptr...) with new sometype(...).
#ifdef __BORLANDC__
		ref_ptr const &operator=(T *other)
		{
			reset(other);
			return *this;
		}
#endif
		///
		/// Get he pointer value, it may return NULL in case of empty pointer
		///
		T *get() const
		{
			return p;
		}
		///
		/// Cast to boolean type: check if the pointer is not empty in similar way as you check ordinary pointers.
		///
		operator bool() const
		{
			return p!=0;
		}
		///
		/// Returns pointer to object, throws cppdb_error if it is NULL
		///
		T *operator->() const
		{
			if(!p)
				throw cppdb_error("cppdb::ref_ptr: attempt to access an empty object");
			return p;
		}
		///
		/// Returns reference to object, throws cppdb_error if it is NULL
		///
		T &operator*() const
		{
			if(!p)
				throw cppdb_error("cppdb::ref_ptr: attempt to access an empty object");
			return *p;
		}
		///
		/// Reset the pointer with new value - old object is dereferenced new is added.
		///
		void reset(T *v=0)
		{
			if(v==p)
				return;
			if(p) {
				if(p->del_ref() == 0) {
					T::dispose(p);
				}
				p=0;
			}
			if(v) {
				v->add_ref();
			}
			p=v;
		}
	private:
		T *p;
	};

	///
	/// \brief This is a class that implements reference counting and designed to be used with ref_ptr
	///
	class ref_counted {
	public:
		///
		/// Create an object with 0 reference count
		///
		ref_counted() : count_(0)
		{
		}
		///
		/// Virtual destructor - for convenience
		///
		virtual ~ref_counted()
		{
		}
		///
		/// Increase reference count
		///
		long add_ref()
		{
			return ++count_;
		}
		///
		/// Get reference count
		///
		long use_count() const
		{
			long val = count_;
			return val;
		}
		///
		/// Decrease reference count
		///
		long del_ref()
		{
			return --count_;
		}
		///
		/// Delete the object
		///
		static void dispose(ref_counted *p)
		{
			delete p;
		}
	private:
		atomic_counter count_;
	};
} // cppdb
#endif
