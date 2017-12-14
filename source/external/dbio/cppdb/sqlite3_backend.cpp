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
#define CPPDB_DRIVER_SOURCE
#ifdef CPPDB_WITH_SQLITE3
# define CPPDB_SOURCE
#endif
#include <sqlite3.h>

#include <cppdb/backend.h>
#include <cppdb/errors.h>
#include <cppdb/utils.h>

#include <sstream>
#include <limits>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifndef __native_client__
// Boost Headers
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>
#endif 

namespace cppdb {
	namespace sqlite3_backend {
		
		class result : public backend::result {
		public:
			result(sqlite3_stmt *st,sqlite3 *conn) : 
				st_(st),
				conn_(conn),
				column_names_prepared_(false),
				cols_(-1)
			{
				cols_=sqlite3_column_count(st_);
			}
			~result() override 
			{
				st_ = nullptr;
			}
			next_row has_next() override
			{
				return next_row_unknown;
			}
			bool next() override 
			{
				int r = sqlite3_step(st_);
				if(r==SQLITE_DONE)
					return false;
				if(r!=SQLITE_ROW) {
					throw cppdb_error(std::string("sqlite3:") + sqlite3_errmsg(conn_));
				}
				return true;
			}
			template<typename T>
			bool do_fetch(int col,T &v)
			{
				if(do_is_null(col))
					return false;
				if(sqlite3_column_type(st_,col)==SQLITE_NULL)
					return false;
				sqlite3_int64 rv = sqlite3_column_int64(st_,col);
				T tmp;
				if(std::numeric_limits<T>::is_signed) {
					tmp=static_cast<T>(rv);
					if(static_cast<sqlite3_int64>(tmp)!=rv)
						throw bad_value_cast();
				}
				else {
					if(rv < 0)
						throw bad_value_cast();
					auto urv = static_cast<unsigned long long>(rv);
					tmp=static_cast<T>(urv);
					if(static_cast<unsigned long long>(tmp)!=urv)
						throw bad_value_cast();
				}
				v=tmp;
				return true;
			}

			bool fetch(int col,boost::uuids::uuid &v) override 
			{
				if(do_is_null(col))
					return false;
				auto const *txt = (char const *)sqlite3_column_text(st_,col);
				memcpy(&v,txt,16);
				return true;
			}
	
			bool fetch(int col,short &v) override 
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,unsigned short &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,int &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,unsigned &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,long &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,unsigned long &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,long long &v) override
			{
				return do_fetch(col,v);
			}
			bool fetch(int col,unsigned long long &v) override
			{
				return do_fetch(col,v);
			}
			template<typename T>
			bool do_real_fetch(int col,T &v)
			{
				if(do_is_null(col))
					return false;
				v=static_cast<T>(sqlite3_column_double(st_,col));
				return true;
			}
			bool fetch(int col,float &v) override 
			{
				return do_real_fetch(col,v);
			}
			bool fetch(int col,double &v) override
			{
				return do_real_fetch(col,v);
			}
			bool fetch(int col,long double &v) override
			{
				return do_real_fetch(col,v);
			}
			bool fetch(int col,std::string &v) override
			{
				if(do_is_null(col))
					return false;
				auto const *txt = (char const *)sqlite3_column_text(st_,col);
				int size = sqlite3_column_bytes(st_,col);
				v.assign(txt,size);
				return true;
			}
			bool fetch(int col,std::ostream &v) override
			{
				if(do_is_null(col))
					return false;
				auto const *txt = (char const *)sqlite3_column_text(st_,col);
				int size = sqlite3_column_bytes(st_,col);
				v.write(txt,size);
				return true;
			}
			bool fetch(int col,std::tm &v) override
			{
				if(do_is_null(col))
					return false;
				v=parse_time((char const *)(sqlite3_column_text(st_,col)));
				return true;
			}
			bool is_null(int col) override
			{
				return do_is_null(col);
			}
			int cols() override 
			{
				return cols_;
			}
			int name_to_column(std::string const &n) override
			{
				if(!column_names_prepared_) {
					for(int i=0;i<cols_;i++) {
						char const *name = sqlite3_column_name(st_,i);
						if(!name) {
							throw std::bad_alloc();
						}
						column_names_[name]=i;
					}
					column_names_prepared_ = true;
				}
				std::map<std::string,int>::const_iterator p=column_names_.find(n);
				if(p==column_names_.end())
					return -1;
				return p->second;
			}
			std::string column_to_name(int col) override
			{
				check(col);
				char const *name = sqlite3_column_name(st_,col);
				if(!name) {
					throw std::bad_alloc();
				}
				return name;
			}
		private:
			bool do_is_null(int col)
			{
				check(col);
				return sqlite3_column_type(st_,col)==SQLITE_NULL;
			}
			void check(int col)
			{
				if(col < 0 || col >= cols_)
					throw invalid_column();
			}
			sqlite3_stmt *st_;
			sqlite3 *conn_;
			std::map<std::string,int> column_names_;
			bool column_names_prepared_;
			int cols_;
		};

		class statement : public backend::statement {
		public:
			void reset() override
			{
				reset_stat();
				sqlite3_clear_bindings(st_);
			}
			void reset_stat()
			{
				if(!reset_) {
					sqlite3_reset(st_);
					reset_=true;
				}
			}
			void bind(int col,boost::uuids::uuid const &v) override 
			{
				#ifndef __native_client__
				reset_stat();
				std::ostringstream ss;
				std::copy(v.begin(), v.end(), std::ostream_iterator<const unsigned char>(ss));
				std::string tmp = ss.str();
				check_bind(sqlite3_bind_blob(st_,col,tmp.c_str(),tmp.size(),SQLITE_TRANSIENT));
			  #endif
			}
			void bind(int col,std::string const &v) override 
			{
				reset_stat();
				check_bind(sqlite3_bind_text(st_,col,v.c_str(),v.size(),SQLITE_STATIC));
			}
			void bind(int col,char const *s) override
			{
				reset_stat();
				check_bind(sqlite3_bind_text(st_,col,s,-1,SQLITE_STATIC));
			}
			void bind(int col,char const *b,char const *e) override 
			{
				reset_stat();
				check_bind(sqlite3_bind_text(st_,col,b,e-b,SQLITE_STATIC));
			}
			void bind(int col,std::tm const &v) override
			{
				reset_stat();
				std::string tmp = cppdb::format_time(v);
				check_bind(sqlite3_bind_text(st_,col,tmp.c_str(),tmp.size(),SQLITE_TRANSIENT));
			}
			void bind(int col,std::istream &v) override 
			{
				reset_stat();
				// TODO Fix me
				std::ostringstream ss;
				ss<<v.rdbuf();
				std::string tmp = ss.str();
				check_bind(sqlite3_bind_text(st_,col,tmp.c_str(),tmp.size(),SQLITE_TRANSIENT));
			}
			void bind(int col,int v) override 
			{
				reset_stat();
				check_bind(sqlite3_bind_int(st_,col,v));
			}
			template<typename IntType>
			void do_bind(int col,IntType value)
			{
				reset_stat();
				int r;
				if(sizeof(value) > sizeof(int) || (long long)(value) > std::numeric_limits<int>::max())
					r = sqlite3_bind_int64(st_,col,static_cast<sqlite3_int64>(value));
				else
					r = sqlite3_bind_int(st_,col,static_cast<int>(value));
				check_bind(r);
			}
			void bind(int col,unsigned v) override 
			{
				do_bind(col,v);
			}
			void bind(int col,long v) override
			{
				do_bind(col,v);
			}
			void bind(int col,unsigned long v) override
			{
				do_bind(col,v);
			}
			void bind(int col,long long v) override
			{
				do_bind(col,v);
			}
			void bind(int col,unsigned long long v) override
			{
				do_bind(col,v);
			}
			void bind(int col,double v) override
			{
				reset_stat();
				check_bind(sqlite3_bind_double(st_,col,v));
			}
			void bind(int col,long double v) override 
			{
				reset_stat();
				check_bind(sqlite3_bind_double(st_,col,static_cast<double>(v)));
			}
			void bind_null(int col) override
			{
				reset_stat();
				check_bind(sqlite3_bind_null(st_,col));
			}
			result *query() override
			{
				reset_stat();
				reset_ = false;
				return new result(st_,conn_);
			}
			long long sequence_last(std::string const &/*name*/) override
			{
				return sqlite3_last_insert_rowid(conn_);
			}
			void exec() override
			{
				reset_stat();
				reset_ = false;
				int r = sqlite3_step(st_);
				if(r!=SQLITE_DONE) {
					if(r==SQLITE_ROW) {
						throw cppdb_error("Using exec with query!");
					}
					else 
						check_bind(r);
				}
			}
			unsigned long long affected() override
			{
				return sqlite3_changes(conn_);
			}
			std::string const &sql_query() override
			{
				return sql_query_;
			}
			statement(std::string const &query,sqlite3 *conn) :
				st_(nullptr),
				conn_(conn),
				reset_(true),
				sql_query_(query)
			{
				if(sqlite3_prepare_v2(conn_,query.c_str(),query.size(),&st_,nullptr)!=SQLITE_OK)
					throw cppdb_error(sqlite3_errmsg(conn_));
			}
			~statement() override
			{
				sqlite3_finalize(st_);
			}

		private:
			void check_bind(int v)
			{
				if(v==SQLITE_RANGE) {
					throw invalid_placeholder(); 
				}
				if(v!=SQLITE_OK) {
					throw cppdb_error(sqlite3_errmsg(conn_));
				}
			}
			sqlite3_stmt *st_;
			sqlite3 *conn_;
			bool reset_;
			std::string sql_query_;
			std::vector<std::string> uuids_;
		};
		class connection : public backend::connection {
		public:
			connection(connection_info const &ci) :
				backend::connection(ci),
				conn_(nullptr)
			{
				std::string dbname=ci.get("db");
				if(dbname.empty()) {
					throw cppdb_error("sqlite3:database file (db propery) not specified");
				}

				std::string mode = ci.get("mode","create");
				int flags = 0;
				if(mode == "create")
					flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;
				else if(mode == "readonly")
					flags = SQLITE_OPEN_READONLY;
				else if(mode == "readwrite")
					flags = SQLITE_OPEN_READWRITE;
				else {
					throw cppdb_error("sqlite3:invalid mode propery, expected "
								" 'create' (default), 'readwrite' or 'readonly' values");
				}

				std::string vfs = ci.get("vfs");
				char const *cvfs = vfs.empty() ? (char const *)nullptr : vfs.c_str();
				
				int busy = ci.get("busy_timeout",-1);
				
				try {
					if(sqlite3_open_v2(dbname.c_str(),&conn_,flags,cvfs)!=SQLITE_OK) {
						if(conn_ == nullptr) {
							throw cppdb_error("sqlite3:failed to create db object");
						}
						
						throw cppdb_error(std::string("sqlite3:Failed to open connection:")
								+ sqlite3_errmsg(conn_));
						
						if(busy!=-1 && sqlite3_busy_timeout(conn_,busy)!=0) {
							throw cppdb_error(std::string("sqlite3:Failed to set timeout:")
									+ sqlite3_errmsg(conn_));
						}
					}
				}
				catch(...) {
					if(conn_) {
						sqlite3_close(conn_);
						conn_ = nullptr;
					}
					throw;
				}
				
			}
			~connection() override 
			{
				sqlite3_close(conn_);
			}
			void begin() override
			{
				fast_exec("begin");	
			}
			void commit() override 
			{
				fast_exec("commit");
			}
			void rollback() override
			{
				fast_exec("rollback");
			}
			statement *prepare_statement(std::string const &q) override
			{
				return new statement(q,conn_);
			}
			statement *create_statement(std::string const &q) override
			{
				return prepare_statement(q);
			}
			std::string escape(std::string const &s) override
			{
				return escape(s.c_str(),s.c_str()+s.size());
			}
			std::string escape(char const *s) override
			{
				return escape(s,s+strlen(s));
			}
			std::string escape(char const *b,char const *e) override
			{
				std::string result;
				result.reserve(e-b);
				for(;b!=e;b++) {
					char c=*b;
					if(c=='\'')
						result+="''";
					else
						result+=c;
				}
				return result;
			}
			std::string driver() override
			{
				return "sqlite3";
			}
			std::string engine() override
			{
				return "sqlite3";
			}
		private:
			void fast_exec(char const *query)
			{
				if(sqlite3_exec(conn_,query,nullptr,nullptr,nullptr)!=SQLITE_OK) {
					throw cppdb_error(std::string("sqlite3:") + sqlite3_errmsg(conn_));
				}
			}

			sqlite3 *conn_;
		};

	} // sqlite3_backend
} // cppdb

extern "C" {
	CPPDB_DRIVER_API cppdb::backend::connection *cppdb_sqlite3_get_connection(cppdb::connection_info const &cs)
	{
		return new cppdb::sqlite3_backend::connection(cs);
	}
}
