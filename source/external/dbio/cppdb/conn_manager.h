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
#ifndef CPPDB_CONN_MANAGER_H
#define CPPDB_CONN_MANAGER_H

#include <cppdb/defs.h>
#include <cppdb/ref_ptr.h>
#include <cppdb/mutex.h>
#include <map>
#include <string>
#include <memory>

namespace cppdb {
	class pool;
	class connection_info;
	namespace backend {
		class connection;
	}

	///
	/// \brief This class is the major gateway to new connections
	///
	/// It handles connection pools and forwards request to the drivers.
	///
	/// This class member functions are thread safe
	///
	class CPPDB_API connections_manager {
		connections_manager();
// Borland erros on hidden destructors in classes without only static methods.
#ifndef __BORLANDC__
		~connections_manager();
#endif
		connections_manager(connections_manager const &);
		void operator = (connections_manager const &);
	public:
		///
		/// Get a singleton instance of the class
		///
		static connections_manager &instance();
		///
		/// Create a new connection using connection string \a cs
		///
		ref_ptr<backend::connection> open(std::string const &cs);
		///
		/// Create a new connection using parsed connection string \a ci
		///
		ref_ptr<backend::connection> open(connection_info const &ci);
		///
		/// Collect all connections that were not used for long time and close them.
		///
		void gc();
	private:
		struct data;
		std::unique_ptr<data> d;

		mutex lock_;
		typedef std::map<std::string,ref_ptr<pool> > connections_type;
		connections_type connections_;
	};
} // cppdb


#endif
