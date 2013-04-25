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
#ifndef CPPDB_MUTEX_H
#define CPPDB_MUTEX_H

#include <cppdb/defs.h>

namespace cppdb {

	///
	/// \brief mutex class, used internally
	///
	class CPPDB_API mutex {
		mutex(mutex const &);
		void operator=(mutex const  &);
	public:
		class guard;
		/// Create mutex
		mutex();
		/// Destroy mutex
		~mutex();
		/// Lock mutex
		void lock();
		/// Unlock mutex
		void unlock();
	private:
		void *mutex_impl_;
	};

	///
	/// \brief scoped guard for mutex
	///
	class mutex::guard {
		guard(guard const &);
		void operator=(guard const &);
	public:
		/// Create scoped lock
		guard(mutex &m) : m_(&m)
		{
			m_->lock();
		}
		/// unlock the mutex
		~guard()
		{
			m_->unlock();
		}
	private:
		mutex *m_;
	};
}
#endif
