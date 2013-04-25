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
#ifdef CPPDB_WITH_PQ
# define CPPDB_SOURCE
#endif 
#include <libpq-fe.h>
#include <libpq/libpq-fs.h>
#include <cppdb/backend.h>
#include <cppdb/errors.h>
#include <cppdb/utils.h>
#include <cppdb/numeric_util.h>
#include <sstream>
#include <vector>
#include <limits>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

#include <iostream>

// Boost Headers
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid_generators.hpp>

namespace cppdb {
	namespace postgresql {
	
		typedef enum {
			lo_type,
			bytea_type
		} blob_type;

		class pqerror : public cppdb_error {
		public:
			pqerror(char const *msg) : cppdb_error(message(msg)) {}
			pqerror(PGresult *r,char const *msg) : cppdb_error(message(msg,r)) {}
			pqerror(PGconn *c,char const *msg) : cppdb_error(message(msg,c)) {}
			
			static std::string message(char const *msg)
			{
				return std::string("cppdb::posgresql: ") + msg;
			}
			static std::string message(char const *msg,PGresult *r)
			{
				std::string m="cppdb::posgresql: ";
				m+=msg;
				m+=": ";
				m+=PQresultErrorMessage(r);
				return m;
			}
			static std::string message(char const *msg,PGconn *c)
			{
				std::string m="cppdb::posgresql: ";
				m+=msg;
				m+=": ";
				m+=PQerrorMessage(c);
				return m;
			}
		};

		class result : public backend::result {
		public:
			result(PGresult *res,PGconn *conn,blob_type b) :
				res_(res),
				conn_(conn),
				rows_(PQntuples(res)),
				cols_(PQnfields(res)),
				current_(-1),
				blob_(b)
			{
				ss_.imbue(std::locale::classic());
			}
			virtual ~result() 
			{
				PQclear(res_);
			}
			virtual next_row has_next()
			{
				if(current_ + 1 < rows_)
					return next_row_exists;
				else
					return last_row_reached; 

			}
			virtual bool next() 
			{
				current_ ++;
				if(current_ < rows_) {
					return true;
				}
				return false;
			}

			template<typename T>
			bool do_fetch(int col,T &v)
			{
				if(do_isnull(col))
					return false;
				std::string tmp(PQgetvalue(res_,current_,col),PQgetlength(res_,current_,col));
				v=parse_number<T>(tmp,ss_);
				return true;
			}
			virtual bool fetch(int col,boost::uuids::uuid &v) 
			{
				std::string s;
				bool result = fetch(col,s);
				boost::uuids::string_generator gen;
				v = gen(s);
				return result;
			}
			virtual bool fetch(int col,short &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,unsigned short &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,int &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,unsigned &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,long &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,unsigned long &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,long long &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,unsigned long long &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,float &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,double &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,long double &v)
			{
				return do_fetch(col,v);
			}
			virtual bool fetch(int col,std::string &v)
			{
				if(do_isnull(col))
					return false;
				v.assign(PQgetvalue(res_,current_,col),PQgetlength(res_,current_,col));
				return true;
			}
			virtual bool fetch(int col,std::ostream &v)
			{
				if(do_isnull(col))
					return false;
				if(blob_ == bytea_type) {
					unsigned char *val=(unsigned char*)PQgetvalue(res_,current_,col);
					size_t len = 0;
					unsigned char *buf=PQunescapeBytea(val,&len);
					if(!buf) {
						throw bad_value_cast();
					}
					try {
						v.write((char *)buf,len);
					}catch(...) {
						PQfreemem(buf);
						throw;
					}
					PQfreemem(buf);
				}
				else { // oid
					Oid id = 0;
					fetch(col,id);
					if(id==0) {
						throw pqerror("fetching large object failed, oid=0");
					}
					int fd = -1;
					try {
						fd = lo_open(conn_,id,INV_READ | INV_WRITE);
						if(fd < 0)
							throw pqerror(conn_,"Failed opening large object for read");
						char buf[4096];
						for(;;) {
							int n=lo_read(conn_,fd,buf,sizeof(buf));
							if(n < 0)
								throw pqerror(conn_,"Failed reading large object");
							if(n>=0)
								v.write(buf,n);
							if(n < int(sizeof(buf)))
								break;
						}
						int r = lo_close(conn_,fd);
						fd = -1;
						if(r < 0)
							throw pqerror(conn_,"error on close of large object");
					}
					catch(...) {
						if(fd != -1)
							lo_close(conn_,fd);
						throw;
					}
				}
				return true;
			}
			virtual bool fetch(int col,std::tm &v)
			{
				if(do_isnull(col))
					return false;
				v=parse_time(PQgetvalue(res_,current_,col));
				return true;
			}
			virtual bool is_null(int col)
			{
				return do_isnull(col);
			}
			virtual int cols() 
			{
				return cols_;
			}
			virtual int name_to_column(std::string const &n) 
			{
				return PQfnumber(res_,n.c_str());
			}
			virtual std::string column_to_name(int pos)
			{
				char const *name = PQfname(res_,pos);
				if(!name)
					return std::string();
				return name;
			}
		private:

			void check(int c)
			{
				if(c < 0 || c>= cols_)
					throw invalid_column();
			}
			bool do_isnull(int col)
			{
				check(col);
				return PQgetisnull(res_,current_,col);
			}
			PGresult *res_;
			PGconn *conn_;
			int rows_;
			int cols_;
			int current_;
			blob_type blob_;
			std::istringstream ss_;
		};

		class statement : public backend::statement {
		public:
			
			typedef enum {
				null_param,
				text_param,
				binary_param
			} param_type;

			statement(PGconn *conn,std::string const &src_query,blob_type b,unsigned long long prepared_id) :
				res_(0),
				conn_(conn),
				orig_query_(src_query),
				params_(0),
				blob_(b)
			{
				fmt_.imbue(std::locale::classic());

				query_.reserve(src_query.size());
				bool inside_string=false;
				for(unsigned i=0;i<src_query.size();i++) {
					char c=src_query[i];
					if(c=='\'') {
						inside_string = !inside_string;
					}
					if(!inside_string && c=='?') {
						query_+='$';
						params_++;
						fmt_<<params_;
						query_+=fmt_.str();
						fmt_.str(std::string());
						fmt_.clear();
					}
					else {
						query_+=c;
					}
				}
				reset();

				if(prepared_id > 0) {

					fmt_.str(std::string());
					fmt_.clear();
					fmt_<<"cppdb_psqlstmt_" << prepared_id;
					prepared_id_ = fmt_.str();
					fmt_.str(std::string());
					fmt_.clear();

					PGresult *r=PQprepare(conn_,prepared_id_.c_str(),query_.c_str(),0,0);
					try {
						if(!r) {
							throw pqerror("Failed to create prepared statement object!");
						}
						if(PQresultStatus(r)!=PGRES_COMMAND_OK)
							throw pqerror(r,"statement preparation failed");
					}
					catch(...) {
						if(r) PQclear(r);
						throw;
					}
					PQclear(r);
				}
			}
			virtual ~statement()
			{
				try {
					if(res_) {
						PQclear(res_);
						res_ = 0;
					}
					if(!prepared_id_.empty()) {
						std::string stmt = "DEALLOCATE " + prepared_id_;
						res_ = PQexec(conn_,stmt.c_str());
						if(res_)  {
							PQclear(res_);
							res_ = 0;
						}
					}
				}
				catch(...) 
				{
				}
			}
			virtual void reset()
			{
				if(res_) {
					PQclear(res_);
					res_ = 0;
				}
				std::vector<std::string> vals(params_);
				std::vector<size_t> lengths(params_,0);
				std::vector<char const *> pvals(params_,0);
				std::vector<param_type> flags(params_,null_param);
				params_values_.swap(vals);
				params_pvalues_.swap(pvals);
				params_plengths_.swap(lengths);
				params_set_.swap(flags);
			}
			virtual void bind(int col,boost::uuids::uuid const &v) 
			{
				std::ostringstream ss;
				std::copy(v.begin(), v.end(), std::ostream_iterator<const unsigned char>(ss));
				params_values_[col-1]=ss.str();
				params_set_[col-1]=binary_param;
			}
			virtual void bind(int col,std::string const &v)
			{
				bind(col,v.c_str(),v.c_str()+v.size());
			}
			virtual void bind(int col,char const *s)
			{
				bind(col,s,s+strlen(s));
			}
			virtual void bind(int col,char const *b,char const *e)
			{
				check(col);
				params_pvalues_[col-1] = b;
				params_plengths_[col-1] = e-b;
				params_set_[col-1]=text_param;
			}
			virtual void bind(int col,std::tm const &v) 
			{
				check(col);
				params_values_[col-1]=cppdb::format_time(v);
				params_set_[col-1]=text_param;
			}
			virtual void bind(int col,std::istream &in)
			{
				check(col);
				if(blob_ == bytea_type) {
					std::ostringstream ss;
					ss << in.rdbuf();
					params_values_[col-1]=ss.str();
					params_set_[col-1]=binary_param;
				}
				else {
					Oid id = 0;
					int fd = -1;
					try {
						id = lo_creat(conn_, INV_READ|INV_WRITE);
						if(id == 0)
							throw pqerror(conn_,"failed to create large object");
						fd = lo_open(conn_,id,INV_READ | INV_WRITE);
						if(fd < 0)
							throw pqerror(conn_,"failed to open large object for writing");
						char buf[4096];
						for(;;) {
							in.read(buf,sizeof(buf));
							int bytes_read = in.gcount();
							if(bytes_read > 0) {
								int n = lo_write(conn_,fd,buf,bytes_read);
								if(n < 0) {
									throw pqerror(conn_,"failed writing to large object");
								}
							}
							if(bytes_read < int(sizeof(buf)))
								break;
						}
						int r = lo_close(conn_,fd);
						fd=-1;
						if(r < 0)
							throw pqerror(conn_,"error closing large object after write");
						bind(col,id);
					}
					catch(...) {
						if(fd<-1)
							lo_close(conn_,fd);
						if(id!=0)
							lo_unlink(conn_,id);
						throw;
					}
				}
			}
			
			template<typename T>
			void do_bind(int col,T v)
			{
				check(col);
				fmt_.str(std::string());
				fmt_.clear();
				if(!std::numeric_limits<T>::is_integer)
					fmt_ << std::setprecision(std::numeric_limits<T>::digits10+1);
				fmt_ << v;
				params_values_[col-1]=fmt_.str();
				params_set_[col-1]=text_param;
				fmt_.str(std::string());
				fmt_.clear();
			}

			virtual void bind(int col,int v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,unsigned v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,long v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,unsigned long v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,long long v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,unsigned long long v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,double v)
			{
				do_bind(col,v);
			}
			virtual void bind(int col,long double v)
			{
				do_bind(col,v);
			}
			virtual void bind_null(int col)
			{
				check(col);
				params_set_[col-1]=null_param;
				std::string tmp;
				params_values_[col-1].swap(tmp);
			}

			void real_query()
			{
				char const * const *pvalues = 0;
				int *plengths = 0;
				int *pformats = 0;
				std::vector<char const *> values;
				std::vector<int> lengths;
				std::vector<int> formats;
				if(params_>0) {
					values.resize(params_,0);
					lengths.resize(params_,0);
					formats.resize(params_,0);
					for(unsigned i=0;i<params_;i++) {
						if(params_set_[i]!=null_param) {
							if(params_pvalues_[i]!=0) {
								values[i]=params_pvalues_[i];
								lengths[i]=params_plengths_[i];
							}
							else {
								values[i]=params_values_[i].c_str();
								lengths[i]=params_values_[i].size();
							}
							if(params_set_[i]==binary_param) {
								formats[i]=1;
							}
						}
					}
					pvalues=&values.front();
					plengths=&lengths.front();
					pformats=&formats.front();
				}
				if(res_) {
					PQclear(res_);
					res_ = 0;
				}
				if(prepared_id_.empty()) {
					res_ = PQexecParams(
						conn_,
						query_.c_str(),
						params_,
						0, // param types
						pvalues,
						plengths,
						pformats, // format - text
						0 // result format - text
						);
				}
				else {
					res_ = PQexecPrepared(
						conn_,
						prepared_id_.c_str(),
						params_,
						pvalues,
						plengths,
						pformats, // format - text
						0 // result format - text
						);
				}
			}

			virtual result *query() 
			{
				real_query();
				switch(PQresultStatus(res_)){
				case PGRES_TUPLES_OK:
					{
						result *ptr = new result(res_,conn_,blob_);
						res_ = 0;
						return ptr;
					}
					break;
				case PGRES_COMMAND_OK:
					throw pqerror("Statement used instread of query");
					break;
				default:
					throw pqerror(res_,"query execution failed ");
				}
			}
			virtual void exec() 
			{
				real_query();
				switch(PQresultStatus(res_)){
				case PGRES_TUPLES_OK:
					throw pqerror("Query used instread of statement");
					break;
				case PGRES_COMMAND_OK:
					break;
				default:
					throw pqerror(res_,"statement execution failed ");
				}

			}
			virtual long long sequence_last(std::string const &sequence)
			{
				PGresult *res = 0;
				long long rowid = 0;
				try {
					char const * const param_ptr = sequence.c_str();
					res = PQexecParams(	conn_,
								"SELECT currval($1)",
								1, // 1 param
								0, // types
								&param_ptr, // param values
								0, // lengths
								0, // formats
								0 // string format
								);
					if(PQresultStatus(res) != PGRES_TUPLES_OK) {
						throw pqerror(res,"failed to fetch last sequence id");
					}
					char const *val = PQgetvalue(res,0,0);
					if(!val || *val==0)
						throw pqerror("Failed to get value for sequecne id");
					fmt_.str(val);
					fmt_.clear();
					fmt_ >> rowid;
					fmt_.str(std::string());
					fmt_.clear();
				}
				catch(...) {
					if(res) PQclear(res);
					throw;
				}
				PQclear(res);
				return rowid;
			}
			virtual unsigned long long affected() 
			{
				if(res_) {
					char const *s=PQcmdTuples(res_);
					if(!s || !*s)
						return 0;
					unsigned long long rows = 0;
					fmt_.str(s);
					fmt_.clear();
					fmt_ >> rows;
					fmt_.str(std::string());
					fmt_.clear();
					return rows;
				}
				return 0;
			}
			virtual std::string const &sql_query()
			{
				return orig_query_;
			}
		private:
			void check(int col)
			{
				if(col < 1 || col > int(params_))
					throw invalid_placeholder();
			}
			PGresult *res_;
			PGconn *conn_;

			std::string query_;
			std::string orig_query_;
			unsigned params_;
			std::vector<std::string> params_values_;
			std::vector<char const *> params_pvalues_;
			std::vector<size_t> params_plengths_;
			std::vector<param_type> params_set_;
			std::string prepared_id_;
			std::stringstream fmt_;
			blob_type blob_;
		};

		class connection : public backend::connection {
		public:
			void do_simple_exec(char const *s)
			{
				PGresult *r=PQexec(conn_,s);
				try {
					
				}
				catch(...) {
					PQclear(r);
					throw;
				}
				PQclear(r);
			}
			virtual void begin()
			{
				do_simple_exec("begin");
			}
			virtual void commit() 
			{
				do_simple_exec("commit");
			}
			virtual void rollback()
			{
				try {
					do_simple_exec("rollback");
				}
				catch(...) {}
			}
			virtual statement *prepare_statement(std::string const &q)
			{
				return new statement(conn_,q,blob_,++prepared_id_);
			}
			virtual statement *create_statement(std::string const &q)
			{
				return new statement(conn_,q,blob_,0);
			}
			std::string do_escape(char const *b,size_t length)
			{
				std::vector<char> buf(2*length+1);
				size_t len = PQescapeStringConn(conn_,&buf.front(),b,length,0);
				return std::string(&buf.front(),len);
			}
			virtual std::string escape(std::string const &s)
			{
				return do_escape(s.c_str(),s.size());
			}
			virtual std::string escape(char const *s)
			{
				return do_escape(s,strlen(s));
			}
			virtual std::string escape(char const *b,char const *e) 
			{
				return do_escape(b,e-b);
			}
			std::string pq_string(connection_info const &ci)
			{
				std::map<std::string,std::string>::const_iterator p;
				std::string pq_str;
				for(p=ci.properties.begin();p!=ci.properties.end();p++) {
					if(p->first.empty() || p->first[0]=='@')
						continue;
					pq_str+=p->first;
					pq_str+="='";
					pq_str+=escape_for_conn(p->second);
					pq_str+="' ";
				}
				return pq_str;
			}
			std::string escape_for_conn(std::string const &v)
			{
				std::string res;
				res.reserve(v.size());
				for(unsigned i=0;i<v.size();i++) {
					if(v[i]=='\\')
						res+="\\\\";
					else if(v[i]=='\'')
						res+="\\\'";
					else
						res+=v[i];
				}
				return res;
			}
			connection(connection_info const &ci) :
				backend::connection(ci),
				conn_(0),
				prepared_id_(0)
			{
				std::string pq=pq_string(ci);
				std::string blob = ci.get("@blob","lo");
				
				if(blob == "bytea")
					blob_ = bytea_type;
				else if(blob == "lo")
					blob_ = lo_type;
				else 
					throw pqerror("@blob property should be either lo or bytea");

				conn_ = 0;
				try {
					conn_ = PQconnectdb(pq.c_str());
					if(!conn_)
						throw pqerror("failed to create connection object");
					if(PQstatus(conn_)!=CONNECTION_OK)
						throw pqerror(conn_,"failed to connect");
				}
				catch(...) {
					if(conn_) {
						PQfinish(conn_);
						conn_ = 0;
					}
					throw;
				}
			}
			virtual ~connection()
			{
				PQfinish(conn_);
			}
			virtual std::string driver()
			{
				return "postgresql";
			}
			virtual std::string engine()
			{
				return "postgresql";
			}
		private:
			PGconn *conn_;
			unsigned long long prepared_id_;
			blob_type blob_;
		};


	} // backend
} // cppdb


extern "C" {
	CPPDB_DRIVER_API cppdb::backend::connection *cppdb_postgresql_get_connection(cppdb::connection_info const &cs)
	{
		return new cppdb::postgresql::connection(cs);
	}
}
