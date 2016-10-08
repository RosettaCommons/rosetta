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
#include <cppdb/pool.h>
#include <cppdb/backend.h>
#include <cppdb/utils.h>
#include <cppdb/driver_manager.h>

#include <stdlib.h>

namespace cppdb {

	struct pool::data {};

	ref_ptr<pool> pool::create(connection_info const &ci)
	{
		ref_ptr<pool> p = new pool(ci);
		return p;
	}
	ref_ptr<pool> pool::create(std::string const &cs)
	{
		connection_info ci(cs);
		ref_ptr<pool> p = new pool(ci);
		return p;
	}

	pool::pool(connection_info const &ci) :
		limit_(0),
		life_time_(0),
		ci_(ci),
		size_(0)
	{
		limit_ = ci_.get("@pool_size",16);
		life_time_ = ci_.get("@pool_max_idle",600);
	}
		
	pool::~pool()
	{
	}

	ref_ptr<backend::connection> pool::open()
	{
		if(limit_ == 0)
			return driver_manager::instance().connect(ci_);

		ref_ptr<backend::connection> p = get();

		if(!p) {
			p=driver_manager::instance().connect(ci_);
		}
		p->set_pool(this);
		return p;
	}

	// this is thread safe member function
	ref_ptr<backend::connection> pool::get()
	{
		if(limit_ == 0)
			return 0;
		ref_ptr<backend::connection> c;
		pool_type garbage;
		std::time_t now = time(0);
		{
			mutex::guard l(lock_);
			// Nothing there should throw so it is safe
			pool_type::iterator p = pool_.begin(),tmp;
			while(p!=pool_.end()) {
				if(p->last_used + life_time_ < now) {
					tmp=p;
					p++;
					garbage.splice(garbage.begin(),pool_,tmp);
					size_ --;
				}
				else {
					// all is sorted by time
					break;
				}
			}
			if(!pool_.empty()) {
				c = pool_.back().conn;
				pool_.pop_back();
				size_ --;
			}
		}
		return c;
	}
	
	// this is thread safe member function
	void pool::put(backend::connection *c_in)
	{
		std::unique_ptr<backend::connection> c(c_in);
		if(limit_ == 0)
			return;
		pool_type garbage;
		std::time_t now = time(0);
		{
			mutex::guard l(lock_);
			// under lock do all very fast
			if(c.get()) {
				pool_.push_back(entry());
				pool_.back().last_used = now;
				pool_.back().conn = c.release();
				size_ ++;
			}
			
			// Nothing there should throw so it is safe
			
			pool_type::iterator p = pool_.begin(),tmp;
			while(p!=pool_.end()) {
				if(p->last_used + life_time_ < now) {
					tmp=p;
					p++;
					garbage.splice(garbage.begin(),pool_,tmp);
					size_ --;
				}
				else {
					// all is sorted by time
					break;
				}
			}
			// can be at most 1 entry bigger then limit
			if(size_ > limit_) {
				garbage.splice(garbage.begin(),pool_,pool_.begin());
				size_--;
			}
		}
	}
	
	void pool::gc()
	{
		put(0);
	}

	void pool::clear()
	{
		pool_type garbage;
		{
			mutex::guard l(lock_);
			garbage.swap(pool_);
			size_ = 0;
		} // destroy outside mutex scope
	}
}



