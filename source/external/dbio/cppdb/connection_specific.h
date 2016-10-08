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
#ifndef CPPDB_CONNECTION_SPECIFIC_H
#define CPPDB_CONNECTION_SPECIFIC_H

#include <cppdb/defs.h>
#include <memory>

namespace cppdb {
	///
	/// \brief Special abstract object that holds a connection specific data
	///
	/// The user is expected to derive its own object from this class
	/// and save them withing the connection
	///
	class CPPDB_API connection_specific_data {
		connection_specific_data(connection_specific_data const &);
		void operator=(connection_specific_data const &);
	public:
		connection_specific_data();
		virtual ~connection_specific_data();

	private:
		struct data;
		std::unique_ptr<data> d;
	};


} // cppdb

#endif
