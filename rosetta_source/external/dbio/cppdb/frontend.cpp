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
#include <cppdb/frontend.h>
#include <cppdb/backend.h>
#include <cppdb/conn_manager.h>

#include <boost/uuid/uuid.hpp>

namespace cppdb {
	struct result::data {};

	result::result() :
		eof_(false),
		fetched_(false),
		current_col_(0)
	{
	}
	result::result(	ref_ptr<backend::result> res,
			ref_ptr<backend::statement> stat,
			ref_ptr<backend::connection> conn)
	: eof_(false),
	  fetched_(false),
	  current_col_(0),
	  res_(res),
	  stat_(stat),
	  conn_(conn)
	{
	}
	result::result(result const &other) :
		eof_(other.eof_),
		fetched_(other.fetched_),
		current_col_(other.current_col_),
		res_(other.res_),
		stat_(other.stat_),
		conn_(other.conn_)
	{
	}

	result const &result::operator=(result const &other)
	{
		eof_ = other.eof_;
		fetched_ = other.fetched_;
		current_col_ = other.current_col_;
		res_ = other.res_;
		stat_ = other.stat_;
		conn_ = other.conn_;
		return *this;
	}

	result::~result()
	{
		clear();
	}

	int result::cols()
	{
		return res_->cols();
	}

	bool result::next()
	{
		if(eof_)
			return false;
		fetched_=true;
		eof_ = res_->next()==false;
		current_col_ = 0;
		return !eof_;
	}
	
	int result::index(std::string const &n)
	{
		int c = res_->name_to_column(n);
		if(c<0)
			throw invalid_column();
		return c;
	}

	std::string result::name(int col)
	{
		if(col < 0 || col>= cols())
			throw invalid_column();
		return res_->column_to_name(col);
	}

	int result::find_column(std::string const &name)
	{
		int c = res_->name_to_column(name);
		if(c < 0)
			return -1;
		return c;
	}

	void result::rewind_column()
	{
		current_col_ = 0;
	}
	
	bool result::empty()
	{
		if(!res_)
			return true;
		return eof_ || !fetched_;
	}

	void result::clear()
	{
		eof_ = true;
		fetched_ = true;
		res_.reset();
		stat_.reset();
		conn_.reset();
	}

	void result::check()
	{
		if(empty())
			throw empty_row_access();
	}

	bool result::is_null(int col)
	{
		return res_->is_null(col);
	}
	bool result::is_null(std::string const &n)
	{
		return is_null(index(n));
	}

	bool result::fetch(int col,boost::uuids::uuid &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,short &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,unsigned short &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,int &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,unsigned &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,long &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,unsigned long &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,long long &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,unsigned long long &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,float &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,double &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,long double &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,std::string &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,std::tm &v) { return res_->fetch(col,v); }
	bool result::fetch(int col,std::ostream &v) { return res_->fetch(col,v); }

    bool result::fetch(std::string const &n,boost::uuids::uuid &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,short &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,unsigned short &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,int &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,unsigned &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,long &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,unsigned long &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,long long &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,unsigned long long &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,float &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,double &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,long double &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,std::string &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,std::tm &v) { return res_->fetch(index(n),v); }
	bool result::fetch(std::string const &n,std::ostream &v) { return res_->fetch(index(n),v); }

    bool result::fetch(boost::uuids::uuid &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(short &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(unsigned short &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(int &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(unsigned &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(long &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(unsigned long &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(long long &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(unsigned long long &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(float &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(double &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(long double &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(std::string &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(std::tm &v) { return res_->fetch(current_col_++,v); }
	bool result::fetch(std::ostream &v) { return res_->fetch(current_col_++,v); }



	struct statement::data {};

	statement::statement() : placeholder_(1) {}
	statement::~statement()
	{
		stat_.reset();
		conn_.reset();
	}

	statement::statement(statement const &other) :
		placeholder_(other.placeholder_),
		stat_(other.stat_),
		conn_(other.conn_)
	{
	}
	statement const &statement::operator=(statement const &other)
	{
		placeholder_ = other.placeholder_;
		stat_=other.stat_;
		conn_=other.conn_;
		return *this;
	}

	statement::statement(ref_ptr<backend::statement> stat,ref_ptr<backend::connection> conn) :
		placeholder_(1),
		stat_(stat),
		conn_(conn)
	{
	}

	void statement::reset()
	{
		placeholder_ = 1;
		stat_->reset();
	}

	statement &statement::operator<<(std::string const &v)
	{
		return bind(v);
	}
	statement &statement::operator<<(char const *s)
	{
		return bind(s);
	}
	
	statement &statement::operator<<(std::tm const &v)
	{
		return bind(v);
	}
	
	statement &statement::operator<<(std::istream &v)
	{
		return bind(v);
	}

	statement &statement::operator<<(void (*manipulator)(statement &st))
	{
		manipulator(*this);
		return *this;
	}
	result statement::operator<<(result (*manipulator)(statement &st))
	{
		return manipulator(*this);
	}
    statement &statement::bind(boost::uuids::uuid v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(int v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(unsigned v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(long v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(unsigned long v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(long long v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(unsigned long long v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(double v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(long double v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}

	statement &statement::bind(std::string const &v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(char const *s)
	{
		stat_->bind(placeholder_++,s);
		return *this;
	}
	statement &statement::bind(char const *b,char const *e)
	{
		stat_->bind(placeholder_++,b,e);
		return *this;
	}
	statement &statement::bind(std::tm const &v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind(std::istream &v)
	{
		stat_->bind(placeholder_++,v);
		return *this;
	}
	statement &statement::bind_null()
	{
		stat_->bind_null(placeholder_++);
		return *this;
	}


    void statement::bind(int col,boost::uuids::uuid v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,std::string const &v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,char const *s)
	{
		stat_->bind(col,s);
	}
	void statement::bind(int col,char const *b,char const *e)
	{
		stat_->bind(col,b,e);
	}
	void statement::bind(int col,std::tm const &v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,std::istream &v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,int v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,unsigned v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,long v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,unsigned long v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,long long v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,unsigned long long v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,double v)
	{
		stat_->bind(col,v);
	}
	void statement::bind(int col,long double v)
	{
		stat_->bind(col,v);
	}
	void statement::bind_null(int col)
	{
		stat_->bind_null(col);
	}

	long long statement::last_insert_id()
	{
		return stat_->sequence_last(std::string());
	}

	long long statement::sequence_last(std::string const &seq)
	{
		return stat_->sequence_last(seq);
	}
	unsigned long long statement::affected()
	{
		return stat_->affected();
	}

	result statement::row()
	{
		ref_ptr<backend::result> backend_res = stat_->query();
		result res(backend_res,stat_,conn_);
		if(res.next()) {
			if(res.res_->has_next() == backend::result::next_row_exists) {
				throw multiple_rows_query();
			}
		}
		return res;
	}
	
	result statement::query()
	{
		ref_ptr<backend::result> res(stat_->query());
		return result(res,stat_,conn_);
	}
	statement::operator result()
	{
		return query();
	}
	void statement::exec() 
	{
		stat_->exec();
	}

	struct session::data {};

	session::session()
	{
	}
	session::session(session const &other) :
		conn_(other.conn_)
	{
	}
	session const &session::operator=(session const &other)
	{
		conn_ = other.conn_;
		return *this;
	}
	session::session(ref_ptr<backend::connection> conn) : 
		conn_(conn) 
	{
	}
	session::~session()
	{
	}
	session::session(connection_info const &ci)
	{
		open(ci);
	}
	session::session(std::string const &cs)
	{
		open(cs);
	}
	
	void session::open(connection_info const &ci)
	{
		conn_ = connections_manager::instance().open(ci);
	}
	void session::open(std::string const &cs)
	{
		conn_ = connections_manager::instance().open(cs);
	}
	void session::close()
	{
		conn_.reset();
	}

	bool session::is_open()
	{
		return conn_;
	}
	
	statement session::prepare(std::string const &query)
	{
		ref_ptr<backend::statement> stat_ptr(conn_->prepare(query));
		statement stat(stat_ptr,conn_);
		return stat;
	}
	
	statement session::create_statement(std::string const &query)
	{
		ref_ptr<backend::statement> stat_ptr(conn_->get_statement(query));
		statement stat(stat_ptr,conn_);
		return stat;
	}
	
	statement session::create_prepared_statement(std::string const &query)
	{
		ref_ptr<backend::statement> stat_ptr(conn_->get_prepared_statement(query));
		statement stat(stat_ptr,conn_);
		return stat;
	}
	
	statement session::create_prepared_uncached_statement(std::string const &query)
	{
		ref_ptr<backend::statement> stat_ptr(conn_->get_prepared_uncached_statement(query));
		statement stat(stat_ptr,conn_);
		return stat;
	}


	statement session::operator<<(std::string const &q)
	{
		return prepare(q);
	}
	statement session::operator<<(char const *s)
	{
		return prepare(s);
	}
	void session::begin()
	{
		conn_->begin();
	}
	void session::commit()
	{
		conn_->commit();
	}
	void session::rollback()
	{
		conn_->rollback();
	}
	std::string session::escape(char const *b,char const *e)
	{
		return conn_->escape(b,e);
	}
	std::string session::escape(char const *s)
	{
		return conn_->escape(s);
	}
	std::string session::escape(std::string const &s)
	{
		return conn_->escape(s);
	}
	std::string session::driver()
	{
		return conn_->driver();
	}
	std::string session::engine()
	{
		return conn_->engine();
	}
	
	struct transaction::data {};

	transaction::transaction(session &s) :
		s_(&s),
		commited_(false)
	{
		s_->begin();
	}
	
	void transaction::commit()
	{
		s_->commit();
		commited_ =true;
	}
	void transaction::rollback()
	{
		if(!commited_)
			s_->rollback();
		commited_=true;
	}
	transaction::~transaction()
	{
		rollback();
	}
	
	void session::clear_cache()
	{
		conn_->clear_cache();
	}




}  // cppdb
