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
#ifndef CPPDB_BACKEND_H
#define CPPDB_BACKEND_H
#include <iosfwd>
#include <ctime>
#include <string>
#include <memory>
#include <map>
#include <typeinfo>
#include <cppdb/defs.h>
#include <cppdb/errors.h>
#include <cppdb/ref_ptr.h>
#include <cppdb/connection_specific.h>

// Borland errors about unknown pool-type without this include.
#ifdef __BORLANDC__
#include <cppdb/pool.h>
#endif

#ifdef PYROSETTA
#include <cppdb/pool.h>
#endif

#include <boost/uuid/uuid.hpp>

namespace cppdb {
	class connection_info;
// Borland needs pool.h, but not this forward declaration.
#ifndef __BORLANDC__
	class pool;
#endif

	///
	/// \brief This namepace includes all classes required to implement a cppdb SQL backend.
	///
	namespace backend {

		///
		/// \brief This class represents query result.
		///
		/// This object is created by statement::query call, backend developer may assume that this object
		/// will stay alive as long as statement that created it exits, i.e. statement would be destroyed after
		/// result.
		///
		class CPPDB_API result : public ref_counted {
		public:

			///
			/// The flag that defines the information about availability of the next row in result
			///
		    enum next_row {
				last_row_reached, ///< No more rows exits, next() would return false
				next_row_exists,  ///< There are more rows, next() would return true
				next_row_unknown  ///< It is unknown, next() may return either true or false
			};

			///
			/// Check if the next row in the result exists. If the DB engine can't perform
			/// this check without loosing data for current row, it should return next_row_unknown.
			///
			virtual next_row has_next() = 0;
			///
			/// Move to next row. Should be called before first access to any of members. If no rows remain
			/// return false, otherwise return true
			///
			virtual bool next() = 0;
			///
			/// Fetch an UUID value for column \a col starting from 0.
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,boost::uuids::uuid &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,short &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,unsigned short &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,int &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,unsigned &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,long &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,unsigned long &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,long long &v) = 0;
			///
			/// Fetch an integer value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to integer or its range is not supported by the integer type.
			///
			virtual bool fetch(int col,unsigned long long &v) = 0;
			///
			/// Fetch a floating point value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to floating point value.
			///
			virtual bool fetch(int col,float &v) = 0;
			///
			/// Fetch a floating point value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to floating point value.
			///
			virtual bool fetch(int col,double &v) = 0;
			///
			/// Fetch a floating point value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, should throw bad_value_cast() if the underlying data
			/// can't be converted to floating point value.
			///
			virtual bool fetch(int col,long double &v) = 0;
			///
			/// Fetch a string value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, any data should be convertible to
			/// text value (as formatting integer, floating point value or date-time as string).
			///
			virtual bool fetch(int col,std::string &v) = 0;
			///
			/// Fetch a BLOB value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid, any data should be convertible to
			/// BLOB value as text (as formatting integer, floating point value or date-time as string).
			///
			virtual bool fetch(int col,std::ostream &v) = 0;
			///
			/// Fetch a date-time value for column \a col starting from 0.
			/// Returns true if ok, returns false if the column value is NULL and the referenced object should remain unchanged
			///
			/// Should throw invalid_column() \a col value is invalid. If the data can't be converted
			/// to date-time it should throw bad_value_cast()
			///
			virtual bool fetch(int col,std::tm &v) = 0;
			///
			/// Check if the column \a col is NULL starting from 0, should throw invalid_column() if the index out of range
			///
			virtual bool is_null(int col) = 0;
			///
			/// Return the number of columns in the result. Should be valid even without calling next() first time.
			///
			virtual int cols() = 0;
			///
			/// Return the number of columns by its name. Return -1 if the name is invalid
			/// Should be able to work even without calling next() first time.
			///
			virtual int name_to_column(std::string const &) = 0;
			///
			/// Return the column name for column index starting from 0.
			/// Should throw invalid_column() if the index out of range
			/// Should be able to work even without calling next() first time.
			///
			virtual std::string column_to_name(int) = 0;

			result();
			virtual ~result();
		private:
			struct data;
			std::unique_ptr<data> d;
		};

		class statements_cache;

		///
		/// \brief This class represents a statement that can be either executed or queried for result
		///
		class CPPDB_API statement : public ref_counted {
		public:
			// Begin of API

			///
			/// Reset the prepared statement to initial state as before the operation. It is
			/// called by front-end each time before new query() or exec() are called.
			///
			virtual void reset() = 0;
			///
			/// Get the query the statement works with. Return it as is, used as key for statement
			/// caching
			///
			virtual std::string const &sql_query() = 0;

			///
			/// Bind a UUID value to column \a col (starting from 1). You may assume
			/// that the reference remains valid until real call of query() or exec()
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,boost::uuids::uuid const &) = 0;
			///
			/// Bind a text value to column \a col (starting from 1). You may assume
			/// that the reference remains valid until real call of query() or exec()
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,std::string const &) = 0;
			///
			/// Bind a text value to column \a col (starting from 1). You may assume
			/// that the reference remains valid until real call of query() or exec()
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,char const *s) = 0;
			///
			/// Bind a text value to column \a col (starting from 1). You may assume
			/// that the reference remains valid until real call of query() or exec()
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,char const *b,char const *e) = 0;
			///
			/// Bind a date-time value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,std::tm const &) = 0;
			///
			/// Bind a BLOB value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,std::istream &) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,int v) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			/// May throw bad_value_cast() if the value out of supported range by the DB.
			///
			virtual void bind(int col,unsigned v) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			/// May throw bad_value_cast() if the value out of supported range by the DB.
			///
			virtual void bind(int col,long v) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			/// May throw bad_value_cast() if the value out of supported range by the DB.
			///
			virtual void bind(int col,unsigned long v) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			/// May throw bad_value_cast() if the value out of supported range by the DB.
			///
			virtual void bind(int col,long long v) = 0;
			///
			/// Bind an integer value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			/// May throw bad_value_cast() if the value out of supported range by the DB.
			///
			virtual void bind(int col,unsigned long long v) = 0;
			///
			/// Bind a floating point value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,double v) = 0;
			///
			/// Bind a floating point value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind(int col,long double v) = 0;
			///
			/// Bind a NULL value to column \a col (starting from 1).
			///
			/// Should throw invalid_placeholder() if the value of col is out of range. May
			/// ignore if it is impossible to know whether the placeholder exists without special
			/// support from back-end.
			///
			virtual void bind_null(int col) = 0;
			///
			/// Fetch the last sequence generated for last inserted row. May use sequence as parameter
			/// if the database uses sequences, should ignore the parameter \a sequence if the last
			/// id is fetched without parameter.
			///
			/// Should be called after exec() for insert statement, otherwise the behavior is undefined.
			///
			/// MUST throw not_supported_by_backend() if such option is not supported by the DB engine.
			///
			virtual long long sequence_last(std::string const &sequence) = 0;
			///
			/// Return the number of affected rows by last statement.
			///
			/// Should be called after exec(), otherwise behavior is undefined.
			///
			virtual unsigned long long affected() = 0;
			///
			/// Return SQL Query result, MAY throw cppdb_error if the statement is not a query
			///
			virtual result *query() = 0;
			///
			/// Execute a statement, MAY throw cppdb_error if the statement returns results.
			///
			virtual void exec() = 0;

			/// \cond INTERNAL
			// Caching support
			static void dispose(statement *selfp);

			void cache(statements_cache *c);
			statement();
			virtual ~statement() ;
			/// \endcond
		private:
			struct data;
			std::unique_ptr<data> d;
			statements_cache *cache_;
		};

		/// \cond INTERNAL
		class CPPDB_API statements_cache {
			statements_cache(statements_cache const &);
			void operator=(statements_cache const &);
		public:
			statements_cache();
			bool active();
			void set_size(size_t n);
			void put(statement *p_in);
			void clear();
			ref_ptr<statement> fetch(std::string const &q);
			~statements_cache();
		private:
			struct data;
			std::unique_ptr<data> d;
		};

		/// \endcond

		class connection;

		///
		/// \brief This class represents a driver that creates connections for
		/// given connection string, custom drivers can be are installed using this
		/// class
		///
		class CPPDB_API driver : public ref_counted {
			driver(driver const &);
			void operator=(driver const &);
		public:
			driver() {}
			virtual ~driver() {}
			///
			/// Return true if the driver in use (i.e. if there is any open connection exist (connection object)
			/// so it can't be removed from the driver
			///
			virtual bool in_use() = 0;
			///
			/// Create a connection object - should be implemented by driver
			///
			virtual connection *open(connection_info const &cs) = 0;
			///
			/// Create a connection object, generally calls open() but may add some information (as registering objects)
			/// and unregistering them
			///
			virtual connection *connect(connection_info const &cs);
		};

		///
		/// \brief This class represents a driver that can be unloaded from the driver_manager.
		///
		class CPPDB_API loadable_driver : public driver {
			loadable_driver(loadable_driver const &);
			void operator=(loadable_driver const &);
		public:
			loadable_driver() {}
			///
			/// Returns true if any of generated connections still exits
			///
			virtual bool in_use();
			virtual ~loadable_driver() {}

			///
			/// Creates a new connection object and keeps track of them for handing (in_use) correctly
			///
			virtual connection *connect(connection_info const &cs);
		};

		extern "C" {
			///
			/// This function type is the function that is generally resolved from the shared objects when loaded
			///
			typedef cppdb::backend::connection *cppdb_backend_connect_function(connection_info const &ci);
		}


		///
		/// \brief Create a static driver using connection function (usable for statically linking drivers).
		///
		class CPPDB_API static_driver : public driver {
		public:
			///
			/// Typedef of the function pointer that is used for creation of connection objects.
			///
			typedef cppdb_backend_connect_function *connect_function_type;

			///
			/// Create a new driver that creates connection using function \a c
			///
			static_driver(connect_function_type c);
			~static_driver();
			///
			/// Always returns true as this driver cannot be unloaded
			///
			bool in_use();
			///
			/// Create new connection - basically calls the function to create the object
			///
			backend::connection *open(connection_info const &ci);
		private:
			connect_function_type connect_;
		};


		///
		/// \brief this class represents connection to database
		///
		class CPPDB_API connection : public ref_counted {
		public:
			///
			/// Create a new object. Connection information \a info is required
			///
			connection(connection_info const &info);
			virtual ~connection();
			/// \cond INTERNAL
			void set_pool(ref_ptr<pool> p);
			ref_ptr<pool> get_pool();
			void set_driver(ref_ptr<loadable_driver> drv);
			static void dispose(connection *c);
			ref_ptr<statement> prepare(std::string const &q);
			ref_ptr<statement> get_prepared_statement(std::string const &q);
			ref_ptr<statement> get_prepared_uncached_statement(std::string const &q);
			ref_ptr<statement> get_statement(std::string const &q);
			/// \endcond

			// API

			///
			/// Start new isolated transaction. Would not be called
			/// withing other transaction on current connection.
			///
			virtual void begin() = 0;
			///
			/// Commit the transaction, you may assume that is called after begin()
			/// was called.
			///
			virtual void commit() = 0;
			///
			/// Rollback the transaction. MUST never throw!!!
			///
			virtual void rollback() = 0;
			///
			/// Create a prepared statement \a q. May throw if preparation had failed.
			/// Should never return null value.
			///
			virtual statement *prepare_statement(std::string const &q) = 0;
			///
			/// Create a (unprepared) statement \a q. May throw if preparation had failed.
			/// Should never return null value.
			///
			virtual statement *create_statement(std::string const &q) = 0;
			///
			/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
			///
			virtual std::string escape(std::string const &) = 0;
			///
			/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
			///
			virtual std::string escape(char const *s) = 0;
			///
			/// Escape a string for inclusion in SQL query. May throw not_supported_by_backend() if not supported by backend.
			///
			virtual std::string escape(char const *b,char const *e) = 0;
			///
			/// Get the name of the driver, for example sqlite3, odbc
			///
			virtual std::string driver() = 0;
			///
			/// Get the name of the SQL Server, for example sqlite3, mssql, oracle, differs from driver() when
			/// the backend supports multiple databases like odbc backend.
			///
			virtual std::string engine() = 0;

			///
			/// Clear statements cache
			///
			void clear_cache();

			///
			/// Check if session specific preparations are done
			///
			/// For new connections always false
			///
			bool once_called() const;

			///
			/// Set once status - true if called flase
			///
			void once_called(bool v);

			///
			/// Get connection specific data pointer of the type \a type , default 0.
			///
			/// The ownership is not changed
			///
			connection_specific_data *connection_specific_get(std::type_info const &type) const;
			///
			/// Release ownership connection specific data pointer of the type \a type
			///
			connection_specific_data *connection_specific_release(std::type_info const &type);
			///
			/// Remove old connection specific data and set new one for a given
			/// type \a type , the ownership on \a p is transferred to connection.
			///
			void connection_specific_reset(std::type_info const &type,connection_specific_data *p = 0);

			///
			/// Check if this back-end can be recycled for reuse in a pool.
			///
			/// If an exception is thrown during operation on DB this flag is reset
			/// to false by the front-end classes result, statement, session.
			///
			/// Default is true
			///
			bool recyclable();

			///
			/// Set recyclable state of the session. If some problem occurs on connection
			/// that prevents its reuse it should be called with false parameter.
			///
			void recyclable(bool value);

		private:

			struct data;
			std::unique_ptr<data> d;
			statements_cache cache_;
			ref_ptr<loadable_driver> driver_;
			ref_ptr<pool> pool_;
			unsigned default_is_prepared_ : 1;
			unsigned once_called_ : 1;
			unsigned recyclable_ : 1;
			unsigned reserverd_ : 29;
		};

	} // backend
} // cppdb

#endif
