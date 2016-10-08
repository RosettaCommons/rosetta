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
#ifndef CPPDB_POOL_H
#define CPPDB_POOL_H

#include <cppdb/defs.h>
#include <cppdb/ref_ptr.h>
#include <cppdb/mutex.h>
#include <cppdb/utils.h>
#include <memory>
#include <list>


namespace cppdb {
	class connection_info;
	namespace backend {
		class connection;
	}
	
	///
	/// \brief Connections pool, allows to handle multiple connections for specific connection string.
	///
	/// Note \ref connections_manager provide more generic interface and hides pools inside it. So you 
	/// generally should use this class only when you prefer to avoid using some global singleton object.
	///
	/// Unlike \ref connections_manager, it uses pool by default unless its size defined as 0.
	///
	/// All this class member functions are thread safe to use from several threads for the same object
	///
	class CPPDB_API pool : public ref_counted {
		pool();
		pool(pool const &);
		void operator=(pool const &);
		pool(connection_info const &ci);
	public:
		/// Create new pool for \a connection_string
		static ref_ptr<pool> create(std::string const &connection_string);
		/// Create new pool for a parsed connection string \a ci
		static ref_ptr<pool> create(connection_info const &ci);
		
		///
		/// Shortcut of cppdb::ref_ptr<cppdb::pool> as cppdb::pool::pointer.
		///
		/// The pointer that is used to handle pool object
		///
		typedef ref_ptr<pool> pointer;

		~pool();

		///
		/// Get a open a connection, it may be fetched either from pool or new one may be created
		///
		ref_ptr<backend::connection> open();
		///
		/// Collect connections that were not used for a long time (close them)
		///
		void gc();

		///
		/// Remove all connections from the pool
		///
		void clear();

		/// \cond INTERNAL
		void put(backend::connection *c_in);
		/// \endcond
	private:
		ref_ptr<backend::connection> get();

		struct data;
		std::unique_ptr<data> d;
		
		struct entry {
			entry() : last_used(0) {}
			ref_ptr<backend::connection> conn;
			std::time_t last_used;
		};

		typedef std::list<entry> pool_type;
		// non-mutable members
		
		size_t limit_;
		int life_time_;
		connection_info ci_;
		
		// mutex protected begin
		mutex lock_;
		size_t size_;
		pool_type pool_;
		// mutex protected end
		
	};
}


#endif
