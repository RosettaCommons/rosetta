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
#include <cppdb/backend.h>
#include <cppdb/utils.h>
#include <cppdb/pool.h>

#include <map>
#include <list>

namespace cppdb {
	namespace backend {
		//result
		struct result::data {};
		result::result() {}
		result::~result() {}
		
		//statement
		struct statement::data {};

		statement::statement() : cache_(0) 
		{
		}
		statement::~statement()
		{
		}
		void statement::cache(statements_cache *c)
		{
			cache_ = c;
		}

		void statement::dispose(statement *p)
		{
			if(!p)
				return;
			statements_cache *cache = p->cache_;
			p->cache_ = 0;
			if(cache) 
				cache->put(p);
			else
				delete p;
		}
		

		//statements cache//////////////

		struct statements_cache::data {

			data() : 
				size(0),
				max_size(0) 
			{
			}

			struct entry;
			typedef std::map<std::string,entry> statements_type;
			typedef std::list<statements_type::iterator> lru_type;
			struct entry {
				ref_ptr<statement> stat;
				lru_type::iterator lru_ptr;
			};
			
			statements_type statements;
			lru_type lru;
			size_t size;
			size_t max_size;


			void insert(ref_ptr<statement> st)
			{
				statements_type::iterator p;
				if((p=statements.find(st->sql_query()))!=statements.end()) {
					p->second.stat = st;
					lru.erase(p->second.lru_ptr);
					lru.push_front(p);
					p->second.lru_ptr = lru.begin();
				}
				else {
					if(size > 0 && size >= max_size) {
						statements.erase(lru.back());
						lru.pop_back();
						size--;
					}
					std::pair<statements_type::iterator,bool> ins = 
						statements.insert(std::make_pair(st->sql_query(),entry()));
					p = ins.first;
					p->second.stat = st;
					lru.push_front(p);
					p->second.lru_ptr = lru.begin();
					size ++;
				}
			}

			ref_ptr<statement> fetch(std::string const &query)
			{
				ref_ptr<statement> st;
				statements_type::iterator p = statements.find(query);
				if(p==statements.end())
					return st;
				st=p->second.stat;
				lru.erase(p->second.lru_ptr);
				statements.erase(p);
				size --;
				return st;
			}

			void clear()
			{
				lru.clear();
				statements.clear();
				size=0;
			}
		}; // data

		statements_cache::statements_cache() 
		{
		}
		void statements_cache::set_size(size_t n)
		{
			if(n!=0 && !active()) {
				d.reset(new data());
				d->max_size = n;
			}
		}
		void statements_cache::put(statement *p_in)
		{
			if(!active()) {
				delete p_in;
			}
			ref_ptr<statement> p(p_in);
			p->reset();
			d->insert(p);
		}
		ref_ptr<statement> statements_cache::fetch(std::string const &q)
		{
			if(!active())
				return 0;
			return d->fetch(q);
		}
		void statements_cache::clear()
		{
			d->clear();
		}
		statements_cache::~statements_cache()
		{
		}

		bool statements_cache::active()
		{
			return d.get()!=0;
		}

		//////////////
		//connection
		//////////////

		struct connection::data {};
		ref_ptr<statement> connection::prepare(std::string const &q) 
		{
			if(default_is_prepared_)
				return get_prepared_statement(q);
			else
				return get_statement(q);
		}
		
		ref_ptr<statement> connection::get_statement(std::string const &q)
		{
			ref_ptr<statement> st = create_statement(q);
			return st;
		}

		ref_ptr<statement> connection::get_prepared_statement(std::string const &q)
		{
			ref_ptr<statement> st;
			if(!cache_.active()) {
				st = prepare_statement(q);
				return st;
			}
			st = cache_.fetch(q);
			if(!st)
				st = prepare_statement(q);
			st->cache(&cache_);
			return st;
		}

		ref_ptr<statement> connection::get_prepared_uncached_statement(std::string const &q)
		{
			ref_ptr<statement> st = prepare_statement(q);
			return st;
		}



		connection::connection(connection_info const &info) :
			pool_(0)
		{
			int cache_size = info.get("@stmt_cache_size",64);
			if(cache_size > 0) {
				cache_.set_size(cache_size);
			}
			std::string def_is_prep = info.get("@use_prepared","on");
			if(def_is_prep == "on")
				default_is_prepared_ = 1;
			else if(def_is_prep == "off") 
				default_is_prepared_ = 0;
			else
				throw cppdb_error("cppdb::backend::connection: @use_prepared should be either 'on' or 'off'");
		}
		connection::~connection()
		{
		}

		void connection::set_pool(ref_ptr<pool> p)
		{
			pool_ = p;
		}
		void connection::set_driver(ref_ptr<loadable_driver> p)
		{
			driver_ = p;
		}
		void connection::clear_cache()
		{
			cache_.clear();
		}
		
		void connection::dispose(connection *c)
		{
			if(!c)
				return;
			ref_ptr<pool> p = c->pool_;
			c->pool_ = 0;
			if(p)
				p->put(c);
			else {
				c->clear_cache();
				// Make sure that driver would not be
				// destoryed destructor of connection exits
				ref_ptr<loadable_driver> driver = c->driver_;
				delete c;
				driver.reset();
			}
		}
		
		connection *driver::connect(connection_info const &cs)
		{
			return open(cs);
		}
		bool loadable_driver::in_use()
		{
			return use_count() > 1;
		}
		connection *loadable_driver::connect(connection_info const &cs)
		{
			connection *c = open(cs);
			c->set_driver(ref_ptr<loadable_driver>(this));
			return c;
		}

		static_driver::static_driver(connect_function_type c) : connect_(c)
		{
		}
		static_driver::~static_driver()
		{
		}
		bool static_driver::in_use()
		{
			return true;
		}
		backend::connection *static_driver::open(connection_info const &ci)
		{
			return connect_(ci);
		}

	} // backend
} // cppdb


