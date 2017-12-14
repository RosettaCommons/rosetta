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
#define CPPDB_SOURCE
#include <cppdb/conn_manager.h>
#include <cppdb/backend.h>
#include <cppdb/pool.h>
#include <cppdb/driver_manager.h>

namespace cppdb {
	struct connections_manager::data{};
	connections_manager::connections_manager() = default;
// Borland erros on hidden destructors in classes without only static methods.
#ifndef __BORLANDC__
	connections_manager::~connections_manager() = default;
#endif

	connections_manager &connections_manager::instance()
	{
		static connections_manager mgr;
		return mgr;
	}

	namespace { 
		struct init { 
			init() 
			{ connections_manager::instance(); }
		} initializer; 
	}

	ref_ptr<backend::connection> connections_manager::open(std::string const &cs)
	{
		ref_ptr<pool> p;
		/// seems we may be using pool
		if(cs.find("@pool_size")!=std::string::npos) {
			mutex::guard l(lock_);
			auto pool_ptr = connections_.find(cs);
			if(pool_ptr!=connections_.end())
				p = pool_ptr->second;
		}

		if(p) {
			return p->open();
		}
		else {
			connection_info ci(cs);
			return open(ci);
		}
	}
	ref_ptr<backend::connection> connections_manager::open(connection_info const &ci)
	{
		if(ci.get("@pool_size",0)==0) {
			return driver_manager::instance().connect(ci);
		}
		ref_ptr<pool> p;
		{
			mutex::guard l(lock_);
			ref_ptr<pool> &ref_p = connections_[ci.connection_string];
			if(!ref_p) {
				ref_p = pool::create(ci);
			}
			p=ref_p;
		}
		return p->open();
	}
	void connections_manager::gc()
	{
		std::vector<ref_ptr<pool> > pools_;
		pools_.reserve(100);
		{
			mutex::guard l(lock_);
			for(auto & connection : connections_) {
				pools_.push_back(connection.second);
			}
		}
		for(auto & pool : pools_) {
			pool->gc();
		}
		pools_.clear();
		{
			mutex::guard l(lock_);
			for(auto p=connections_.begin();p!=connections_.end();) {
				if(p->second->use_count() == 1) {
					pools_.push_back(p->second);
					auto tmp = p;
					++p;
					connections_.erase(tmp);
				}
				else
					++p;
			}
		}
		pools_.clear();

	}
		

} // cppdb
