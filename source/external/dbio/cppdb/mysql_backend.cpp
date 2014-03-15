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
#ifdef CPPDB_WITH_MYSQL
# define CPPDB_SOURCE
#endif
#include <mysql/mysql.h>

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

namespace mysql_backend {

class cppdb_myerror : public cppdb_error
{
public:
	cppdb_myerror(std::string const &str) :
		cppdb_error("cppdb::mysql::" + str)
	{
	}
};
namespace unprep {
	class result : public backend::result {
	public:

		///
		/// Check if the next row in the result exists. If the DB engine can't perform
		/// this check without loosing data for current row, it should return next_row_unknown.
		///
		virtual next_row has_next()
		{
			if(!res_)
				return last_row_reached;
			if(current_row_ >= mysql_num_rows(res_))
				return last_row_reached;
			else
				return next_row_exists;
		}
		///
		/// Move to next row. Should be called before first access to any of members. If no rows remain
		/// return false, otherwise return true
		///
		virtual bool next()
		{
			if(!res_)
				return false;
			current_row_ ++;
			row_ = mysql_fetch_row(res_);
			if(!row_)
				return false;
			return true;
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		char const *at(int col)
		{
			if(!res_)
				throw empty_row_access();
			if(col < 0 || col >= cols_)
				throw invalid_column();
			return row_[col];
		}

		char const *at(int col,size_t &len)
		{
			if(!res_)
				throw empty_row_access();
			if(col < 0 || col >= cols_)
				throw invalid_column();
			unsigned long *lengths = mysql_fetch_lengths(res_);
			if(lengths==0)
				throw cppdb_myerror("Can't get length of column");
			len = lengths[col];
			return row_[col];
		}

		template<typename T>
		bool do_fetch(int col,T &v)
		{
			size_t len;
			char const *s=at(col,len);
			if(!s)
				return false;
			v = parse_number<T>(std::string(s,len),fmt_);
			return true;
		}
		virtual bool fetch(int col,boost::uuids::uuid &v) {
			size_t len;
			char const *s=at(col,len);
			if(!s)
				return false;
			memcpy(&v,s,len);
			return true;
		}
		virtual bool fetch(int col,short &v)
		{
			return do_fetch(col,v);;
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
		///
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
		///
		virtual bool fetch(int col,std::string &v)
		{
			size_t len;
			char const *s=at(col,len);
			if(!s)
				return false;
			v.assign(s,len);
			return true;
		}
		virtual bool fetch(int col,std::ostream &v)
		{
			size_t len;
			char const *s=at(col,len);
			if(!s)
				return false;
			v.write(s,len);
			return true;
		}
		///
		/// Fetch a date-time value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid. If the data can't be converted
		/// to date-time it should throw bad_value_cast()
		///
		virtual bool fetch(int col,std::tm &v)
		{
			size_t len;
			char const *s=at(col,len);
			if(!s)
				return false;
			std::string tmp(s,len);
			v = parse_time(tmp);
			return true;
		}
		///
		/// Check if the column \a col is NULL starting from 0, should throw invalid_column() if the index out of range
		///
		virtual bool is_null(int col)
		{
			return at(col) == 0;
		}
		///
		/// Return the number of columns in the result. Should be valid even without calling next() first time.
		///
		virtual int cols()
		{
			return cols_;
		}
		virtual std::string column_to_name(int col)
		{
			if(col < 0 || col >=cols_)
				throw invalid_column();
			if(!res_)
				throw empty_row_access();
			MYSQL_FIELD *flds=mysql_fetch_fields(res_);
			if(!flds) {
				throw cppdb_myerror("Internal error empty fileds");
			}
			return flds[col].name;
		}
		virtual int name_to_column(std::string const &name)
		{
			if(!res_)
				throw empty_row_access();
			MYSQL_FIELD *flds=mysql_fetch_fields(res_);
			if(!flds) {
				throw cppdb_myerror("Internal error empty fileds");
			}
			for(int i=0;i<cols_;i++)
				if(name == flds[i].name)
					return i;
			return -1;
		}

		// End of API

		result(MYSQL *conn) :
			res_(0),
			cols_(0),
			current_row_(0),
			row_(0)
		{
			fmt_.imbue(std::locale::classic());
			res_ = mysql_store_result(conn);
			if(!res_) {
				cols_ = mysql_field_count(conn);
				if(cols_ == 0)
					throw cppdb_myerror("Seems that the query does not produce any result");
			}
			else {
				cols_ = mysql_num_fields(res_);
			}

		}
		~result()
		{
			if(res_)
				mysql_free_result(res_);
		}
	private:
		std::istringstream fmt_;
		MYSQL_RES *res_;
		int cols_;
		unsigned current_row_;
		MYSQL_ROW row_;
	};

	class statement : public backend::statement {
	public:
		virtual std::string const &sql_query()
		{
			return query_;
		}

		virtual void bind(int col,boost::uuids::uuid const &v)
		{
			std::ostringstream ss;
			std::copy(v.begin(), v.end(), std::ostream_iterator<const unsigned char>(ss));
			std::string tmp=ss.str();
			bind(col,tmp);
		}

		virtual void bind(int col,std::string const &s)
		{
			bind(col,s.c_str(),s.c_str()+s.size());
		}
		virtual void bind(int col,char const *s)
		{
			bind(col,s,s+strlen(s));
		}
		virtual void bind(int col,char const *b,char const *e)
		{
			std::vector<char> buf(2*(e-b)+1);
			size_t len = mysql_real_escape_string(conn_,&buf.front(),b,e-b);
			std::string &s=at(col);
			s.clear();
			s.reserve(e-b+2);
			s+='\'';
			s.append(&buf.front(),len);
			s+='\'';
		}
		virtual void bind(int col,std::tm const &v)
		{
			std::string &s = at(col);
			s.clear();
			s.reserve(30);
			s+='\'';
			s+=cppdb::format_time(v);
			s+='\'';
		}
		virtual void bind(int col,std::istream &v)
		{
			std::ostringstream ss;
			ss << v.rdbuf();
			std::string tmp=ss.str();
			bind(col,tmp);
		}
		template<typename T>
		void do_bind(int col,T v)
		{
			fmt_.str(std::string());
			if(!std::numeric_limits<T>::is_integer)
				fmt_ << std::setprecision(std::numeric_limits<T>::digits10+1);
			fmt_ << v;
			std::string tmp = fmt_.str();
			at(col).swap(tmp);
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
			at(col)="NULL";
		}
		virtual long long sequence_last(std::string const &/*sequence*/)
		{
			return mysql_insert_id(conn_);
		}
		virtual unsigned long long affected()
		{
			return mysql_affected_rows(conn_);
		}

		void bind_all(std::string &real_query)
		{
			size_t total = query_.size();
			for(unsigned i=0;i<params_.size();i++) {
				total+=params_[i].size();
			}
			real_query.clear();
			real_query.reserve(total);
			size_t pos_ = 0;
			for(unsigned i=0;i<params_.size();i++) {
				size_t marker = binders_[i];
				real_query.append(query_,pos_,marker-pos_);
				pos_ = marker+1;
				real_query.append(params_[i]);
			}
			real_query.append(query_,pos_,std::string::npos);
		}

		virtual result *query()
		{
			std::string real_query;
			bind_all(real_query);
			reset_params();
			if(mysql_real_query(conn_,real_query.c_str(),real_query.size())) {
				throw cppdb_myerror(mysql_error(conn_));
			}
			return new result(conn_);
		}

		virtual void exec()
		{
			std::string real_query;
			bind_all(real_query);
			reset_params();
			if(mysql_real_query(conn_,real_query.c_str(),real_query.size())) {
				throw cppdb_myerror(mysql_error(conn_));
			}
			MYSQL_RES *r=mysql_store_result(conn_);
			if(r){
				mysql_free_result(r);
				throw cppdb_myerror("Calling exec() on query!");
			}
		}

		std::string &at(int col)
		{
			if(col < 1 || col > params_no_) {
				throw invalid_placeholder();
			}
			return params_[col-1];
		}
		void reset_params()
		{
			params_.clear();
			params_.resize(params_no_,"NULL");
		}

		statement(std::string const &q,MYSQL *conn) :
			query_(q),
			conn_(conn),
			params_no_(0)
		{
			fmt_.imbue(std::locale::classic());
			bool inside_text = false;
			for(size_t i=0;i<query_.size();i++) {
				if(query_[i]=='\'') {
					inside_text=!inside_text;
				}
				if(query_[i]=='?' && !inside_text) {
					params_no_++;
					binders_.push_back(i);
				}
			}
			if(inside_text) {
				throw cppdb_myerror("Unterminated string found in query");
			}
			reset_params();
		}
		virtual ~statement()
		{
		}
		virtual void reset()
		{
		}

	private:
		std::ostringstream fmt_;
		std::vector<std::string> params_;
		std::vector<size_t> binders_;

		std::string query_;
		MYSQL *conn_;
		int params_no_;
	};
} // uprep

namespace prep {

	class result : public backend::result {
		struct bind_data {
			bind_data() :
				ptr(0),
				length(0),
				is_null(0),
				error(0)
			{
				memset(&buf,0,sizeof(buf));
			}
			char buf[128];
			std::vector<char> vbuf;
			char *ptr;
			unsigned long length;
			my_bool is_null;
			my_bool error;
		};
	public:

		///
		/// Check if the next row in the result exists. If the DB engine can't perform
		/// this check without loosing data for current row, it should return next_row_unknown.
		///
		virtual next_row has_next()
		{
			if(current_row_ >= mysql_stmt_num_rows(stmt_))
				return last_row_reached;
			else
				return next_row_exists;
		}
		///
		/// Move to next row. Should be called before first access to any of members. If no rows remain
		/// return false, otherwise return true
		///
		virtual bool next()
		{
			current_row_ ++;
			reset();
			if(cols_ > 0) {
				if(mysql_stmt_bind_result(stmt_,&bind_[0])) {
					throw cppdb_myerror(mysql_stmt_error(stmt_));
				}
			}
			int r = mysql_stmt_fetch(stmt_);
			if(r==MYSQL_NO_DATA) {
				return false;
			}
			if(r==MYSQL_DATA_TRUNCATED) {
				for(int i=0;i<cols_;i++) {
					if(bind_data_[i].error && !bind_data_[i].is_null && bind_data_[i].length >= sizeof(bind_data_[i].buf)) {
						bind_data_[i].vbuf.resize(bind_data_[i].length);
						bind_[i].buffer = &bind_data_[i].vbuf.front();
						bind_[i].buffer_length = bind_data_[i].length;
						if(mysql_stmt_fetch_column(stmt_,&bind_[i],i,0)) {
							throw cppdb_myerror(mysql_stmt_error(stmt_));
						}
						bind_data_[i].ptr = &bind_data_[i].vbuf.front();
					}
				}
			}
			return true;
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		bind_data &at(int col)
		{
			if(col < 0 || col >= cols_)
				throw invalid_column();
			if(bind_data_.empty())
				throw cppdb_myerror("Attempt to access data without fetching it first");
			return bind_data_.at(col);
		}

		template<typename T>
		bool do_fetch(int col,T &v)
		{
			bind_data &d=at(col);
			if(d.is_null)
				return false;

			v=parse_number<T>(std::string(d.ptr,d.length),fmt_);
			return true;
		}
		virtual bool fetch(int col,boost::uuids::uuid &v)
		{
			bind_data &d=at(col);
			if(d.is_null)
				return false;
			memcpy(&v,d.ptr,d.length);
			return true;
		}
		virtual bool fetch(int col,short &v)
		{
			return do_fetch(col,v);;
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,unsigned short &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,int &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,unsigned &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,long &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,unsigned long &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,long long &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch an integer value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to integer or its range is not supported by the integer type.
		///
		virtual bool fetch(int col,unsigned long long &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch a floating point value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to floating point value.
		///
		virtual bool fetch(int col,float &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch a floating point value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to floating point value.
		///
		virtual bool fetch(int col,double &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch a floating point value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
		/// can't be converted to floating point value.
		///
		virtual bool fetch(int col,long double &v)
		{
			return do_fetch(col,v);
		}
		///
		/// Fetch a string value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, any data should be convertible to
		/// text value (as formatting integer, floating point value or date-time as string).
		///
		virtual bool fetch(int col,std::string &v)
		{
			bind_data &d=at(col);
			if(d.is_null)
				return false;
			v.assign(d.ptr,d.length);
			return true;
		}
		///
		/// Fetch a BLOB value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid, any data should be convertible to
		/// BLOB value as text (as formatting integer, floating point value or date-time as string).
		///
		virtual bool fetch(int col,std::ostream &v)
		{
			bind_data &d=at(col);
			if(d.is_null)
				return false;
			v.write(d.ptr,d.length);
			return true;
		}
		///
		/// Fetch a date-time value for column \a col starting from 0.
		/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
		///
		/// Should throw invalid_column() \a col value is invalid. If the data can't be converted
		/// to date-time it should throw bad_value_cast()
		///
		virtual bool fetch(int col,std::tm &v)
		{
			std::string tmp;
			if(!fetch(col,tmp))
				return false;
			v = parse_time(tmp);
			return true;
		}
		///
		/// Check if the column \a col is NULL starting from 0, should throw invalid_column() if the index out of range
		///
		virtual bool is_null(int col)
		{
			return at(col).is_null;
		}
		///
		/// Return the number of columns in the result. Should be valid even without calling next() first time.
		///
		virtual int cols()
		{
			return cols_;
		}
		virtual std::string column_to_name(int col)
		{
			if(col < 0 || col >=cols_)
				throw invalid_column();
			MYSQL_FIELD *flds=mysql_fetch_fields(meta_);
			if(!flds) {
				throw cppdb_myerror("Internal error empty fileds");
			}
			return flds[col].name;
		}
		virtual int name_to_column(std::string const &name)
		{
			MYSQL_FIELD *flds=mysql_fetch_fields(meta_);
			if(!flds) {
				throw cppdb_myerror("Internal error empty fileds");
			}
			for(int i=0;i<cols_;i++)
				if(name == flds[i].name)
					return i;
			return -1;
		}

		// End of API

		result(MYSQL_STMT *stmt) :
			stmt_(stmt), current_row_(0),meta_(0)
		{
			fmt_.imbue(std::locale::classic());
			cols_ = mysql_stmt_field_count(stmt_);
			if(mysql_stmt_store_result(stmt_)) {
				throw cppdb_myerror(mysql_stmt_error(stmt_));
			}
			meta_ = mysql_stmt_result_metadata(stmt_);
			if(!meta_) {
				throw cppdb_myerror("Seems that the query does not produce any result");
			}
		}
		~result()
		{
			mysql_free_result(meta_);
		}
		void reset()
		{
			bind_.resize(0);
			bind_data_.resize(0);
			bind_.resize(cols_,MYSQL_BIND());
			bind_data_.resize(cols_,bind_data());
			for(int i=0;i<cols_;i++) {
				bind_[i].buffer_type = MYSQL_TYPE_STRING;
				bind_[i].buffer = bind_data_[i].buf;
				bind_[i].buffer_length = sizeof(bind_data_[i].buf);
				bind_[i].length = &bind_data_[i].length;
				bind_[i].is_null = &bind_data_[i].is_null;
				bind_[i].error = &bind_data_[i].error;
				bind_data_[i].ptr = bind_data_[i].buf;
			}
		}
	private:
		std::istringstream fmt_;
		int cols_;
		MYSQL_STMT *stmt_;
		unsigned current_row_;
		MYSQL_RES *meta_;
		std::vector<MYSQL_BIND> bind_;
		std::vector<bind_data> bind_data_;
	};

	class statement : public backend::statement {
		struct param {
			my_bool is_null;
			bool is_blob;
			unsigned long length;
			std::string value;
			void *buffer;
			param() :
				is_null(1),
				is_blob(false),
				length(0),
				buffer(0)
			{
			}
			void set(char const *b,char const *e,bool blob=false)
			{
				length = e - b;
				buffer = const_cast<char *>(b);
				is_blob = blob;
				is_null = 0;
			}
			void set_str(std::string const &s)
			{
				value = s;
				buffer = const_cast<char *>(value.c_str());
				length = value.size();
				is_null = 0;
			}
			void set(std::tm const &t)
			{
				set_str(cppdb::format_time(t));
			}
			void bind_it(MYSQL_BIND *b)
			{
				b->is_null = &is_null;
				if(!is_null) {
					b->buffer_type = is_blob ? MYSQL_TYPE_BLOB : MYSQL_TYPE_STRING;
					b->buffer = buffer;
					b->buffer_length = length;
					b->length = &length;
				}
				else {
					b->buffer_type = MYSQL_TYPE_NULL;
				}
			}
		};

	public:
		// Begin of API
		///
		/// Get the query the statement works with. Return it as is, used as key for statement
		/// caching
		///
		virtual std::string const &sql_query()
		{
			return query_;
		}

		///
		/// Bind a UUID value to column \a col (starting from 1). You may assume
		/// that the reference remains valid until real call of query() or exec()
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
    virtual void bind(int col,boost::uuids::uuid const &v)
		{
			std::ostringstream ss;
			std::copy(v.begin(), v.end(), std::ostream_iterator<const unsigned char>(ss));
			at(col).set_str(ss.str());
			at(col).is_blob = true;
		}

		///
		/// Bind a text value to column \a col (starting from 1). You may assume
		/// that the reference remains valid until real call of query() or exec()
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,std::string const &s)
		{
			at(col).set(s.c_str(),s.c_str()+s.size());
		}
		///
		/// Bind a text value to column \a col (starting from 1). You may assume
		/// that the reference remains valid until real call of query() or exec()
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,char const *s)
		{
			at(col).set(s,s+strlen(s));
		}
		///
		/// Bind a text value to column \a col (starting from 1). You may assume
		/// that the reference remains valid until real call of query() or exec()
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,char const *b,char const *e)
		{
			at(col).set(b,e);
		}
		///
		/// Bind a date-time value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,std::tm const &v)
		{
			at(col).set(v);
		}
		///
		/// Bind a BLOB value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,std::istream &v)
		{
			std::ostringstream ss;
			ss << v.rdbuf();
			at(col).set_str(ss.str());
			at(col).is_blob = true;
		}
		template<typename T>
		void do_bind(int col,T v)
		{
			fmt_.str(std::string());
			if(!std::numeric_limits<T>::is_integer)
				fmt_ << std::setprecision(std::numeric_limits<T>::digits10+1);
			fmt_ << v;
			at(col).set_str(fmt_.str());
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,int v)
		{
			do_bind(col,v);
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		/// May throw bad_value_cast() if the value out of supported range by the DB.
		///
		virtual void bind(int col,unsigned v)
		{
			do_bind(col,v);
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		/// May throw bad_value_cast() if the value out of supported range by the DB.
		///
		virtual void bind(int col,long v)
		{
			do_bind(col,v);
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		/// May throw bad_value_cast() if the value out of supported range by the DB.
		///
		virtual void bind(int col,unsigned long v)
		{
			do_bind(col,v);
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		/// May throw bad_value_cast() if the value out of supported range by the DB.
		///
		virtual void bind(int col,long long v)
		{
			do_bind(col,v);
		}
		///
		/// Bind an integer value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		/// May throw bad_value_cast() if the value out of supported range by the DB.
		///
		virtual void bind(int col,unsigned long long v)
		{
			do_bind(col,v);
		}
		///
		/// Bind a floating point value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,double v)
		{
			do_bind(col,v);
		}
		///
		/// Bind a floating point value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind(int col,long double v)
		{
			do_bind(col,v);
		}
		///
		/// Bind a NULL value to column \a col (starting from 1).
		///
		/// Should throw invalid_placeholder() if the value of col is out of range. May
		/// ignore if it is impossible to know whether the placeholder exists without special
		/// support from back-end.
		///
		virtual void bind_null(int col)
		{
			at(col)=param();
		}
		///
		/// Fetch the last sequence generated for last inserted row. May use sequence as parameter
		/// if the database uses sequences, should ignore the parameter \a sequence if the last
		/// id is fetched without parameter.
		///
		/// Should be called after exec() for insert statement, otherwise the behavior is undefined.
		///
		/// MUST throw not_supported_by_backend() if such option is not supported by the DB engine.
		///
		virtual long long sequence_last(std::string const &/*sequence*/)
		{
			return mysql_stmt_insert_id(stmt_);
		}
		///
		/// Return the number of affected rows by last statement.
		///
		/// Should be called after exec(), otherwise behavior is undefined.
		///
		virtual unsigned long long affected()
		{
			return mysql_stmt_affected_rows(stmt_);
		}

		void bind_all()
		{
			if(!params_.empty()) {
				for(unsigned i=0;i<params_.size();i++)
					params_[i].bind_it(&bind_[i]);
				if(mysql_stmt_bind_param(stmt_,&bind_.front())) {
					throw cppdb_myerror(mysql_stmt_error(stmt_));
				}
			}
		}

		///
		/// Return SQL Query result, MAY throw cppdb_error if the statement is not a query
		///
		virtual result *query()
		{
			bind_all();
			if(mysql_stmt_execute(stmt_)) {
				throw cppdb_myerror(mysql_stmt_error(stmt_));
			}
			return new result(stmt_);
		}
		///
		/// Execute a statement, MAY throw cppdb_error if the statement returns results.
		///
		virtual void exec()
		{
			bind_all();
			if(mysql_stmt_execute(stmt_)) {
				if(strstr(mysql_stmt_error(stmt_),"Deadlock"))
				{
					throw backend_deadlock();
				}else
				{
					throw cppdb_myerror(mysql_stmt_error(stmt_));
				}
			}
			if(mysql_stmt_store_result(stmt_)) {
				if(strstr(mysql_stmt_error(stmt_),"Deadlock"))
				{
					throw backend_deadlock();
				}else
				{
					throw cppdb_myerror(mysql_stmt_error(stmt_));
				}
			}
			MYSQL_RES *r=mysql_stmt_result_metadata(stmt_);
			if(r) {
				mysql_free_result(r);
				throw cppdb_myerror("Calling exec() on query!");
			}
		}
		// End of API

		// Caching support

		statement(std::string const &q,MYSQL *conn) :
			query_(q),
			stmt_(0),
			params_count_(0)
		{
			fmt_.imbue(std::locale::classic());

			stmt_ = mysql_stmt_init(conn);
			try {
				if(!stmt_) {
					throw cppdb_myerror(" Failed to create a statement");
				}
				if(mysql_stmt_prepare(stmt_,q.c_str(),q.size())) {
					throw cppdb_myerror(mysql_stmt_error(stmt_));
				}
				params_count_ = mysql_stmt_param_count(stmt_);
				reset_data();
			}
			catch(...) {
				if(stmt_)
					mysql_stmt_close(stmt_);
				throw;
			}
		}
		virtual ~statement()
		{
			mysql_stmt_close(stmt_);
		}
		void reset_data()
		{
			params_.resize(0);
			params_.resize(params_count_);
			bind_.resize(0);
			bind_.resize(params_count_,MYSQL_BIND());
		}
		virtual void reset()
		{
			reset_data();
			mysql_stmt_reset(stmt_);
		}

	private:
		param &at(int col)
		{
			if(col < 1 || col > params_count_)
				throw invalid_placeholder();
			return params_[col-1];
		}

		std::ostringstream fmt_;
		std::vector<param> params_;
		std::vector<MYSQL_BIND> bind_;
		std::string query_;
		MYSQL_STMT *stmt_;
		int params_count_;
	};

} // prep

class connection;

class connection : public backend::connection {
public:
	connection(connection_info const &ci) :
		backend::connection(ci),
		conn_(0)
	{
		conn_ = mysql_init(0);
		if(!conn_) {
			throw cppdb_error("cppdb::mysql failed to create connection");
		}
		std::string host = ci.get("host","");
		char const *phost = host.empty() ? 0 : host.c_str();
		std::string user = ci.get("user","");
		char const *puser = user.empty() ? 0 : user.c_str();
		std::string password = ci.get("password","");
		char const *ppassword = password.empty() ? 0 : password.c_str();
		std::string database = ci.get("database","");
		char const *pdatabase = database.empty() ? 0 : database.c_str();
		int port = ci.get("port",0);
		std::string unix_socket = ci.get("unix_socket","");
		char const *punix_socket = unix_socket.empty() ? 0 : unix_socket.c_str();

#if MYSQL_VERSION_ID >= 50507
		std::string default_auth = ci.get("default_auth","");
		if (!default_auth.empty()) {
			mysql_set_option(MYSQL_DEFAULT_AUTH, default_auth.c_str());
		}
#endif
		std::string init_command = ci.get("init_command","");
		if(!init_command.empty()) {
			mysql_set_option(MYSQL_INIT_COMMAND, init_command.c_str());
		}
		if(ci.has("opt_compress")) {
			if(ci.get("opt_compress", 1)) {
				mysql_set_option(MYSQL_OPT_COMPRESS, NULL);
			}
		}
		if(ci.has("opt_connect_timeout")) {
			if(unsigned connect_timeout = ci.get("opt_connect_timeout", 0)) {
				mysql_set_option(MYSQL_OPT_CONNECT_TIMEOUT, &connect_timeout);
			}
		}
		if(ci.has("opt_guess_connection")) {
			if(ci.get("opt_guess_connection", 1)) {
				mysql_set_option(MYSQL_OPT_GUESS_CONNECTION, NULL);
			}
		}
		if(ci.has("opt_local_infile")) {
			if(unsigned local_infile = ci.get("opt_local_infile", 0)) {
				mysql_set_option(MYSQL_OPT_CONNECT_TIMEOUT, &local_infile);
			}
		}
		if(ci.has("opt_named_pipe")) {
			if(ci.get("opt_named_pipe", 1)) {
				mysql_set_option(MYSQL_OPT_NAMED_PIPE, NULL);
			}
		}
		if(ci.has("opt_protocol")) {
			if(unsigned protocol = ci.get("opt_protocol", 0)) {
				mysql_set_option(MYSQL_OPT_PROTOCOL, &protocol);
			}
		}
		if(ci.has("opt_read_timeout")) {
			if(unsigned read_timeout = ci.get("opt_read_timeout", 0)) {
				mysql_set_option(MYSQL_OPT_READ_TIMEOUT, &read_timeout);
			}
		}
		if(ci.has("opt_reconnect")) {
			if(unsigned reconnect = ci.get("opt_reconnect", 1)) {
				my_bool value = reconnect;
				mysql_set_option(MYSQL_OPT_RECONNECT, &value);
			}
		}
#if MYSQL_VERSION_ID >= 50507
		std::string plugin_dir = ci.get("plugin_dir", "");
		if(!plugin_dir.empty()) {
			mysql_set_option(MYSQL_PLUGIN_DIR, plugin_dir.c_str());
		}
#endif
		std::string set_client_ip = ci.get("set_client_ip", "");
		if(!set_client_ip.empty()) {
			mysql_set_option(MYSQL_SET_CLIENT_IP, set_client_ip.c_str());
		}
		if(ci.has("opt_ssl_verify_server_cert")) {
			if(unsigned verify = ci.get("opt_ssl_verify_server_cert", 1)) {
				my_bool value = verify;
				mysql_set_option(MYSQL_OPT_SSL_VERIFY_SERVER_CERT, &value);
			}
		}
		if(ci.has("opt_use_embedded_connection")) {
			if(ci.get("opt_use_embedded_connection", 1)) {
				mysql_set_option(MYSQL_OPT_USE_EMBEDDED_CONNECTION, NULL);
			}
		}
		if(ci.has("opt_use_remote_connection")) {
			if(ci.get("opt_use_remote_connection", 1)) {
				mysql_set_option(MYSQL_OPT_USE_REMOTE_CONNECTION, NULL);
			}
		}
		if(ci.has("opt_write_timeout")) {
			if(unsigned write_timeout = ci.get("opt_write_timeout", 0)) {
				mysql_set_option(MYSQL_OPT_WRITE_TIMEOUT, &write_timeout);
			}
		}
		std::string read_default_file = ci.get("read_default_file", "");
		if(!read_default_file.empty()) {
			mysql_set_option(MYSQL_READ_DEFAULT_FILE, read_default_file.c_str());
		}
		std::string read_default_group = ci.get("read_default_group", "");
		if(!read_default_group.empty()) {
			mysql_set_option(MYSQL_READ_DEFAULT_GROUP, read_default_group.c_str());
		}
		if(ci.has("report_data_truncation")) {
			if(unsigned report = ci.get("report_data_truncation", 1)) {
				my_bool value = report;
				mysql_set_option(MYSQL_REPORT_DATA_TRUNCATION, &value);
			}
		}
#if MYSQL_VERSION_ID >= 40101
		if(ci.has("secure_auth")) {
			if(unsigned secure = ci.get("secure_auth", 1)) {
				my_bool value = secure;
				mysql_set_option(MYSQL_SECURE_AUTH, &value);
			}
		}
#endif
		std::string set_charset_dir = ci.get("set_charset_dir", "");
		if(!set_charset_dir.empty()) {
			mysql_set_option(MYSQL_SET_CHARSET_DIR, set_charset_dir.c_str());
		}
		std::string set_charset_name = ci.get("set_charset_name", "");
		if(!set_charset_name.empty()) {
			mysql_set_option(MYSQL_SET_CHARSET_NAME, set_charset_name.c_str());
		}
		std::string shared_memory_base_name = ci.get("shared_memory_base_name", "");
		if(!shared_memory_base_name.empty()) {
			mysql_set_option(MYSQL_SHARED_MEMORY_BASE_NAME, shared_memory_base_name.c_str());
		}

		if(!mysql_real_connect(conn_,phost,puser,ppassword,pdatabase,port,punix_socket,0)) {
			std::string err="unknown";
			try { err = mysql_error(conn_); }catch(...){}
			mysql_close(conn_);
			throw cppdb_myerror(err);
		}

	}
	~connection()
	{
		mysql_close(conn_);
	}
	// API

	void exec(std::string const &s)
	{
		if(mysql_real_query(conn_,s.c_str(),s.size())) {
			throw cppdb_myerror(mysql_error(conn_));
		}
	}

	///
	/// Start new isolated transaction. Would not be called
	/// withing other transaction on current connection.
	///
	virtual void begin()
	{
		exec("BEGIN");
	}
	///
	/// Commit the transaction, you may assume that is called after begin()
	/// was called.
	///
	virtual void commit()
	{
		exec("COMMIT");
	}
	///
	/// Rollback the transaction. MUST never throw!!!
	///
	virtual void rollback()
	{
		try {
			exec("ROLLBACK");
		}
		catch(...) {
		}
	}
	///
	/// Create a prepared statement \a q. May throw if preparation had failed.
	/// Should never return null value.
	///
	virtual backend::statement *prepare_statement(std::string const &q)
	{
		return new prep::statement(q,conn_);
	}
	virtual backend::statement *create_statement(std::string const &q)
	{
		return new unprep::statement(q,conn_);
	}
	///
	/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
	///
	virtual std::string escape(std::string const &s)
	{
		return escape(s.c_str(),s.c_str()+s.size());
	}
	///
	/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
	///
	virtual std::string escape(char const *s)
	{
		return escape(s,s+strlen(s));
	}
	///
	/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
	///
	virtual std::string escape(char const *b,char const *e)
	{
		std::vector<char> buf(2*(e-b)+1);
		size_t len = mysql_real_escape_string(conn_,&buf.front(),b,e-b);
		std::string result;
		result.assign(&buf.front(),len);
		return result;
	}
	///
	/// Get the name of the driver, for example sqlite3, odbc
	///
	virtual std::string driver()
	{
		return "mysql";
	}
	///
	/// Get the name of the SQL Server, for example sqlite3, mssql, oracle, differs from driver() when
	/// the backend supports multiple databases like odbc backend.
	///
	virtual std::string engine()
	{
		return "mysql";
	}

	// API

private:
	///
	/// Set a custom MYSQL option on the connection.
	///
	///
	void mysql_set_option(mysql_option option, const void* arg)
	{
		// char can be casted to void but not the other way, support older API
		if(mysql_options(conn_, option, static_cast<char const *>(arg))) {
			throw cppdb_error("cppdb::mysql failed to set option");
		}
	}
	connection_info ci_;
	MYSQL *conn_;
};


} // backend
} // cppdb

extern "C" {
	CPPDB_DRIVER_API cppdb::backend::connection *cppdb_mysql_get_connection(cppdb::connection_info const &cs)
	{
		return new cppdb::mysql_backend::connection(cs);
	}
}
