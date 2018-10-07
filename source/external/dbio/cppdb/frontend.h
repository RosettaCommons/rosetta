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
#ifndef CPPDB_FRONTEND_H
#define CPPDB_FRONTEND_H
#include <cppdb/defs.h>
#include <cppdb/errors.h>
#include <cppdb/ref_ptr.h>

// Borland errors about unknown pool-type without this include.
#ifdef __BORLANDC__
#include <cppdb/backend.h>
#endif

#ifdef PYROSETTA
#include <cppdb/backend.h>
#endif

#include <iosfwd>
#include <ctime>
#include <string>
#include <memory>
#include <typeinfo>

#include <boost/uuid/uuid.hpp>

///
/// The namespace of all data related to the cppdb api
///

namespace cppdb {

	class result;
	class statement;
	class session;
	class connection_info;
	class connection_specific_data;

	///
	/// Get CppDB Version String. It consists of "A.B.C", where A
	/// is a major version number, B is a minor version number and
	/// C is patch version
	///
	CPPDB_API char const *version_string();
	///
	/// Return CppDB version as a number as a sum A * 10000 + B * 100 + C
	/// where A is a major version number, B is a minor version number and
	/// C is patch version
	///
	CPPDB_API int version_number();

// Borland needs pool.h, but not this forward declaration.
#ifndef __BORLANDC__
	namespace backend {
		class result;
		class statement;
		class connection;
	}
#endif

	///
	/// Null value marker
	///
	typedef enum {
		null_value,  	///< The value is null value
		not_null_value	///< The valus is not a null value
	} null_tag_type;

	/// \cond INTERNAL
	namespace tags {
		template<typename T>
		struct into_tag {
			T &value;
			null_tag_type &tag;
			into_tag(T &v, null_tag_type &t) : value(v),tag(t) {}
		};

		template<typename T>
		struct use_tag {
			T value;
			null_tag_type tag;
			use_tag(T v,null_tag_type t) : value(v),tag(t) {}
		};

	} // tags
	/// \endcond

	///
	/// \brief  Create a pair of value and tag for fetching a value from row.
	///
	/// The fetched
	/// value will be stored in \a value if the column is not null and the flag
	/// if the value is null or not saved in \a tag
	///
	template<typename T>
	tags::into_tag<T> into(T &value,null_tag_type &tag)
	{
		return tags::into_tag<T>(value,tag);
	}

	///
	/// \brief Create a pair of a string value and tag for storing it to DB
	///

	inline tags::use_tag<std::string const &> use(std::string const &v,null_tag_type tag)
	{
		return tags::use_tag<std::string const &>(v,tag);
	}

	///
	/// \brief Create a pair of a string value and tag for storing it to DB
	///

	inline tags::use_tag<char const *> use(char const *v,null_tag_type tag)
	{
		return tags::use_tag<char const *>(v,tag);
	}

	///
	/// \brief Create a pair of value and tag for storing it to DB
	///
	template<typename T>
	tags::use_tag<T> use(T value,null_tag_type tag)
	{
		return tags::use_tag<T>(value,tag);
	}

	/// \cond INTERNAL
	namespace details {
		template<typename Object>
		class functor {
		public:
			functor(functor const &other)  :
				functor_(other.functor_),
				wrapper_(other.wrapper_)
			{
			}
			functor const &operator=(functor const &other)
			{
				functor_ = other.functor_;
				wrapper_ = other.wrapper_;
				return *this;
			}
			functor(void (*func)(Object &))
			{
				functor_ = reinterpret_cast<void *>(reinterpret_cast<size_t>(func));
				wrapper_ = &functor::call_func;

			}
			template<typename RealFunctor>
			functor(RealFunctor const &f)
			{
				// The usual casts are not enough for all compilers
				functor_ = reinterpret_cast<void const *>(&f);
				wrapper_ = &functor<Object>::template call_it<RealFunctor>;
			}
			void operator()(Object &p) const
			{
				wrapper_(functor_,p);
			}
		private:
			static void call_func(void const *pointer,Object &parameter)
			{
				typedef void function_type(Object &);
				function_type *f = reinterpret_cast<function_type *>(reinterpret_cast<size_t>((pointer)));
				f(parameter);
			}
			template<typename Functor>
			static void call_it(void const *pointer,Object &parameter)
			{
				Functor const *f_ptr = reinterpret_cast<Functor const *>(pointer);
				Functor const &f=*f_ptr;
				f(parameter);
			}
			void const *functor_;
			void (*wrapper_)(void const *,Object &);
		};
	} // details
	/// \endcond

	#ifdef CPPDB_DOXYGEN
	///
	/// Special object that can be constructed from generic function like object \a f.
	///
	/// So once_functor(f) can be created if f can be used like f(s) where s is \ref cppdb::session
	///
	typedef unspecified_class_type once_functor;
	#else
	typedef details::functor<session> once_functor;
	#endif

	///
	/// \brief This object represents query result.
	///
	/// This object and it is generally created by statement::query() call, default constructor
	/// is provided for consistency, but access to any member function with exception of empty() would
	/// throw an exception.
	///
	class CPPDB_API result {
	public:
		///
		/// Create an empty result, it is not useful except for having default constructor
		///
		result();
		///
		/// Destroys the result, note, if the result of statement is not destroyed, it would
		/// not be returned to statements cache.
		///
		~result();
		///
		/// Copy result, note it only keeps the reference to actual object so copy is just
		/// copy of the reference
		///
		result(result const &);
		///
		/// Assign result, note it only keeps the reference to actual object so assignment is just
		/// copy of the reference
		///
		result const &operator=(result const &);

		///
		/// Return the number of columns in the result
		///
		int cols();
		///
		/// Move forward to next row, returns false if no more rows available.
		///
		/// Notes:
		///
		/// - You should call next() at least once before you use fetch() functions
		/// - You must not call fetch() functions if next() returned false, it would cause empty_row_access exception.
		///
		bool next();

		///
		/// Convert column name \a n to its index, throws invalid_column if the name is not valid.
		///
		int index(std::string const &n);
		///
		/// Convert column name \a n to its index, returns -1 if the name is not valid.
		///
		int find_column(std::string const &name);

		///
		/// Convert column index to column name, throws invalid_column if col is not in range 0<= col < cols()
		///
		std::string name(int col);

		///
		/// Return true if the column number \a col (starting from 0) has NULL value
		///
		bool is_null(int col);
		///
		/// Return true if the column named \a n has NULL value
		///
		bool is_null(std::string const &n);

		///
		/// Clears the result, no further use of the result should be done until it is assigned again with a new statement result.
		///
		/// It is useful when you want to release all data and return the statement to cache
		///
		void clear();
		///
		/// Reset current column index, so fetch without column index can be used once again
		///
		void rewind_column();
		///
		/// Check if the current row is empty, it is in 3 cases:
		///
		/// -# Empty result
		/// -# next() wasn't called first time
		/// -# next() returned false;
		///
		bool empty();


		///
		/// Fetch a value from column \a col (starting from 0) into \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same it tries to cast the data, if casting fails or the
		/// data is out of the type range, throws bad_value_cast().
		///
		bool fetch(int col,boost::uuids::uuid &v);
		///
		/// \copydoc fetch(int,uuid&)
		///
		bool fetch(int col,short &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,unsigned short &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,int &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,unsigned &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,long &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,unsigned long &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,long long &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,unsigned long long &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,float &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,double &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,long double &v);
		///
		/// Fetch a textual value from column \a col (starting from 0) into \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same, if possible it converts it into textual representation.
		///
		bool fetch(int col,std::string &v);
		///
		/// \copydoc fetch(int,short&)
		///
		bool fetch(int col,std::tm &v);
		///
		/// Fetch a binary large object value from column \a col (starting from 0) into a stream \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not blob, it may throw bad_value_cast()
		///
		bool fetch(int col,std::ostream &v);

		///
		/// Fetch a value from column named \a n into \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same it tries to cast the data, if casting fails or the
		/// data is out of the type range, throws bad_value_cast().
		///
		/// If the \a n value is invalid throws invalid_column exception
		///
		bool fetch(std::string const &n,boost::uuids::uuid &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,short &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,unsigned short &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,int &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,unsigned &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,long &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,unsigned long &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,long long &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,unsigned long long &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,float &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,double &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,long double &v);
		///
		/// Fetch a textual value from column named \a n into \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same, if possible it converts it into textual representation. If
		/// the \a n value is invalid throws invalid_column exception
		///
		bool fetch(std::string const &n,std::string &v);
		///
		/// \copydoc fetch(std::string const &,short&)
		///
		bool fetch(std::string const &n,std::tm &v);
		///
		/// Fetch a binary large object value from column named \a name into a stream \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not blob, it may throw bad_value_cast(). If
		/// the \a n value is invalid throws invalid_column exception
		///
		bool fetch(std::string const &n,std::ostream &v);


		///
		/// Fetch a value from the next column in the row starting from the first one. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same it tries to cast the data, if casting fails or the
		/// data is out of the type range, throws bad_value_cast().
		///
		/// If fetch was called more times then cols() it throws invalid_column exception, to use
		/// it once again from the beginning on the same row call rewind_column() member function.
		/// It is not required to call rewind_column() after calling next() as column index is reset
		/// automatically.
		///
		bool fetch(boost::uuids::uuid &v);
		/// \copydoc fetch(uuid&)
		bool fetch(short &v);
		/// \copydoc fetch(short&)
		bool fetch(unsigned short &v);
		/// \copydoc fetch(short&)
		bool fetch(int &v);
		/// \copydoc fetch(short&)
		bool fetch(unsigned &v);
		/// \copydoc fetch(short&)
		bool fetch(long &v);
		/// \copydoc fetch(short&)
		bool fetch(unsigned long &v);
		/// \copydoc fetch(short&)
		bool fetch(long long &v);
		/// \copydoc fetch(short&)
		bool fetch(unsigned long long &v);
		/// \copydoc fetch(short&)
		bool fetch(float &v);
		/// \copydoc fetch(short&)
		bool fetch(double &v);
		/// \copydoc fetch(short&)
		bool fetch(long double &v);
		///
		/// Fetch a textual value from the next column in the row starting from the first one. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not same, if possible it converts it into textual representation.
		///
		/// If fetch was called more times then cols() it throws invalid_column exception, to use
		/// it once again from the beginning on the same row call rewind_column() member function.
		/// It is not required to call rewind_column() after calling next() as column index is reset
		/// automatically.
		///
		bool fetch(std::string &v);
		/// \copydoc fetch(short&)
		bool fetch(std::tm &v);
		///
		/// Fetch a blob value from the next column in the row starting from the first one into stream \a v. Returns false
		/// if the value in NULL and \a v is not updated, otherwise returns true.
		///
		/// If the data type is not blob, it may throw bad_value_cast().
		///
		/// If fetch was called more times then cols() it throws invalid_column exception, to use
		/// it once again from the beginning on the same row call rewind_column() member function.
		/// It is not required to call rewind_column() after calling next() as column index is reset
		/// automatically.
		///
		bool fetch(std::ostream &v);

		///
		/// Get a value of type \a T from column named \a name. If the column
		/// is null throws null_value_fetch(), if the column \a name is invalid throws invalid_column,
		/// if the column value cannot be converted to type T (see fetch functions) it throws bad_value_cast.
		///

		template<typename T>
		T get(std::string const &name)
		{
			T v=T();
			if(!fetch(name,v))
				throw null_value_fetch();
			return v;
		}

		///
		/// Get a value of type \a T from column named \a name. If the column
		/// is null returns \a def, if the column \a name is invalid throws invalid_column,
		/// if the column value cannot be converted to type T (see fetch functions) it throws bad_value_cast.
		///
		template<typename T>
		T get(std::string const &name, T const &def)
		{
			T v=T();
			if(!fetch(name,v))
				return def;
			return v;
		}

		///
		/// Get a value of type \a T from column \a col (starting from 0). If the column
		/// is null throws null_value_fetch(), if the column index is invalid throws invalid_column,
		/// if the column value cannot be converted to type T (see fetch functions) it throws bad_value_cast.
		///
		template<typename T>
		T get(int col)
		{
			T v=T();
			if(!fetch(col,v))
				throw null_value_fetch();
			return v;
		}

		///
		/// Get a value of type \a T from column \a col (starting from 0). If the column
		/// is null returns \a def, if the column index is invalid throws invalid_column,
		/// if the column value cannot be converted to type T (see fetch functions) it throws bad_value_cast.
		///
		template<typename T>
		T get(int col, T const &def)
		{
			T v=T();
			if(!fetch(col,v))
				return def;
			return v;
		}

		///
		/// Syntactic sugar, used together with into() function.
		///
		/// res << into(x,y) is same as
		///
		/// \code
		/// y = res.fetch(x) ? not_null_value : null_value
		/// \endcode
		///
		template<typename T>
		result &operator>>(tags::into_tag<T> ref)
		{
			if(fetch(ref.value))
				ref.tag = not_null_value;
			else
				ref.tag = null_value;
			return *this;
		}

		///
		/// Syntactic sugar, same as fetch(\a value)
		///
		template<typename T>
		result &operator>>(T &value)
		{
			fetch(value);
			return *this;
		}


	private:
		result(	ref_ptr<backend::result> res,
			ref_ptr<backend::statement> stat,
			ref_ptr<backend::connection> conn);

		void check();

		friend class statement;

		struct data;

		std::unique_ptr<data> d;

		bool eof_;
		bool fetched_;
		int current_col_;
		ref_ptr<backend::result> res_;
		ref_ptr<backend::statement> stat_;
		ref_ptr<backend::connection> conn_;
	};

	///
	/// \brief This class represents a prepared (or ordinary) statement that can be executed.
	///
	/// This object is usually created via session::prepare() function.
	///
	class CPPDB_API statement {
	public:
		///
		/// Default constructor, provided for convenience, access to any member function
		/// of empty statement will cause an exception being thrown.
		///
		statement();
		///
		/// Destructor, it releases prepared statement and if the statements cache is enabled
		/// it returns it into the cache.
		///
		/// Note: if result object created by this statement is alive, the underlying backned::statement would
		/// be on hold until result object is destroyed and only then the statement would be put back
		/// into cache.
		///
		~statement();
		///
		/// Copy statement.
		///
		/// Please note it copies only the reference to underlying statement object, so copies
		/// of same statement represent same object and it is strongly not recommended to access the underlying
		/// backend::statement by two different statement objects.
		///
		statement(statement const &);
		///
		/// Assign a statement.
		///
		/// Please note it copies only the reference to underlying statement object, so copies
		/// of same statement represent same object and it is strongly not recommended to access the underlying
		/// backend::statement by two different statement objects.
		///
		statement const &operator=(statement const &);

		///
		/// Reset the statement - remove all bindings and return it into initial state so query() or exec()
		/// functions can be called once again.
		///
		/// You must use it if you use the same statement multiple times.
		///
		/// Note, it is different from clear() where the statement is fully released and access to
		/// it would throw an exception
		///
		void reset();

		///
		/// Clear the statement, removes it, any access to statement object would throw an exception till it would
		/// be assigned once again
		///

		void clear();

		///
		/// Check if the statement is empty, it is empty when created with default constructor or when cleared
		/// with clear() member function.
		///
		bool empty() const;

		///
		/// Bind a value \a v to the next placeholder (starting from the first) marked with '?' marker in the query.
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind(boost::uuids::uuid v);
		/// \copydoc bind(uuid)
		statement &bind(int v);
		/// \copydoc bind(int)
		statement &bind(unsigned v);
		/// \copydoc bind(int)
		statement &bind(long v);
		/// \copydoc bind(int)
		statement &bind(unsigned long v);
		/// \copydoc bind(int)
		statement &bind(long long v);
		/// \copydoc bind(int)
		statement &bind(unsigned long long v);
		/// \copydoc bind(int)
		statement &bind(double v);
		/// \copydoc bind(int)
		statement &bind(long double v);
		///
		/// Bind a string value \a v to the next placeholder marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind(std::string const &v);
		///
		/// Bind a null terminated string value \a s to the next placeholder marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind(char const *s);
		///
		/// Bind a string value in range [\a b, \a e ) to the next placeholder marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind(char const *b,char const *e);
		/// \copydoc bind(int)
		statement &bind(std::tm const &v);
		///
		/// Bind a BLOB value \a v to the next placeholder marked with '?' marker in the query.
		///
		/// Note: the reference to the stream MUST remain valid until the statement is queried or executed!
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind(std::istream &v);
		///
		/// Bind a NULL value  to the next placeholder marked with '?' marker in the query.
		///
		/// If number of calls is higher then the number placeholders is the statement it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		statement &bind_null();

// Without the following statement &operator<<(T v) errors for tags::use_tag<T> as T.
#ifdef __BORLANDC__
		template<typename T>
		statement &bind(tags::use_tag<T> const &val)
		{
			if(val.tag == null_value)
				return bind_null();
			else
				return bind(val.value);
		}
#endif

		///
		/// Bind a value \a v to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind(int col,boost::uuids::uuid v);
		/// \copydoc bind(int,uuid)
		void bind(int col,int v);
		/// \copydoc bind(int,int)
		void bind(int col,unsigned v);
		/// \copydoc bind(int,int)
		void bind(int col,long v);
		/// \copydoc bind(int,int)
		void bind(int col,unsigned long v);
		/// \copydoc bind(int,int)
		void bind(int col,long long v);
		/// \copydoc bind(int,int)
		void bind(int col,unsigned long long v);
		/// \copydoc bind(int,int)
		void bind(int col,double v);
		/// \copydoc bind(int,int)
		void bind(int col,long double v);
		///
		/// Bind a string value \a v to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind(int col,std::string const &v);
		///
		/// Bind a null terminated string value \a s to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind(int col,char const *s);
		///
		/// Bind a string value in range [\a b, \a e ) to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// Note: the reference to the string MUST remain valid until the statement is queried or executed!
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind(int col,char const *b,char const *e);
		/// \copydoc bind(int,int)
		void bind(int col,std::tm const &v);
		///
		/// Bind a BLOB value \a v to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// Note: the reference to the stream MUST remain valid until the statement is queried or executed!
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind(int col,std::istream &v);
		///
		/// Bind a NULL value to the placeholder number \a col (starting from 1) marked with '?' marker in the query.
		///
		/// If \a cols is invalid (less then 1 or higher then the number of the placeholders is the statement) it
		/// may throw invalid_placeholder exception.
		///
		/// If placeholder was not binded the behavior is undefined and may vary between different backends.
		///
		void bind_null(int col);

		///
		/// Get last insert id from the last executed statement, note, it is the same as sequence_last("").
		///
		/// Some backends requires explicit sequence name so you should use sequence_last("sequence_name") in such
		/// case.
		///
		/// If the statement is actually query, the behavior is undefined and may vary between backends.
		///
		long long last_insert_id();
		///
		/// Get last created sequence value from the last executed statement.
		///
		/// If the backend does not support named sequences but rather supports "auto increment" columns (like MySQL, Sqlite3),
		/// the \a seq parameter is ignored.
		///
		///
		/// If the statement is actually query, the behavior is undefined and may vary between backends.
		///
		long long sequence_last(std::string const &seq);
		///
		/// Get the number of affected rows by the last statement,
		///
		///
		/// If the statement is actually query, the behavior is undefined and may vary between backends.
		///
		unsigned long long affected();

		///
		/// Fetch a single row from the query. Unlike query(), you should not call result::next()
		/// function as it is already called. You may check if the data was fetched using result::empty()
		/// function.
		///
		/// If the result set consists of more then one row it throws multiple_rows_query exception, however some backends
		/// may ignore this.
		///
		/// If the statement is not query statement (like SELECT) it would likely
		/// throw an exception, however the behavior may vary between backends that may ignore this error.
		///
		result row();
		///
		/// Fetch a result of the query, if the statement is not query statement (like SELECT) it would likely
		/// throw an exception, however the behavior may vary between backends that may ignore this error.
		///
		result query();
		///
		/// Same as query() - syntactic sugar
		///
		operator result();

		///
		/// Execute a statement, of the statement is actually SELECT like operator, it throws cppdb_error exception,
		/// however the behavior may vary between backends that may ignore this error.
		///
		void exec();

		///
		/// Same as bind(v);
		///
		statement &operator<<(std::string const &v);
		///
		/// Same as bind(s);
		///
		statement &operator<<(char const *s);
		///
		/// Same as bind(v);
		///
		statement &operator<<(std::tm const &v);
		///
		/// Same as bind(v);
		///
		statement &operator<<(std::istream &v);
		///
		/// Apply manipulator on the statement, same as manipulator(*this).
		///
		statement &operator<<(void (*manipulator)(statement &st));
		///
		/// Apply manipulator on the statement, same as manipulator(*this).
		///
		result operator<<(result (*manipulator)(statement &st));

		///
		/// Used together with use() function.
		///
		/// The call st<<use(x,tag) is same as
		///
		/// \code
		///  (tag == null_value) ?  st.bind_null() : st.bind(x)
		/// \endcode
		///
		template<typename T>
		statement &operator<<(tags::use_tag<T> const &val)
		{
			if(val.tag == null_value)
				return bind_null();
			else
				return bind(val.value);
		}

		///
		/// Same as bind(v);
		///
		template<typename T>
		statement &operator<<(T v)
		{
			return bind(v);
		}

	private:
		statement(ref_ptr<backend::statement> stat,ref_ptr<backend::connection> conn);

		friend class session;

		int placeholder_;
		ref_ptr<backend::statement> stat_;
		ref_ptr<backend::connection> conn_;
		struct data;
		std::unique_ptr<data> d;
	};

	///
	/// \brief Manipulator that causes statement execution. Used as:
	///
	/// \code
	///  sql << "delete from test" << cppdb::exec;
	/// \endcode
	///
	inline void exec(statement &st)
	{
		st.exec();
	}

	///
	/// \brief Manipulator that binds null value. Used as:
	///
	/// \code
	///  sql << "insert into foo values(?,?,?)" << x << cppdb::null << y << cppdb::exec;
	/// \endcode
	///
	inline void null(statement &st)
	{
		st.bind_null();
	}

	///
	/// \brief Manipulator that fetches a single row. Used as:
	///
	/// \code
	///  cppdb::result r = sql << "SELECT name where uid=?" << id << cppdb::row;
	/// if(!r.empty()) {
	///  ...
	/// }
	/// \endcode
	///
	/// Or:
	///
	/// \code
	///  sql << "SELECT name where uid=?" << id << cppdb::row >> name;
	/// \endcode
	///
	/// Which would throw empty_row_access exception on attempt to fetch name if the result is empty.
	///
	inline result row(statement &st)
	{
		return st.row();
	}

	///
	/// \brief SQL session object that represents a single connection and is the gateway to SQL database
	///
	/// It is the main class that is used for access to the DB, it uses various singleton classes to
	/// load drivers open connections and cache them.
	///
	class CPPDB_API session {
	public:

		///
		/// Create an empty session object, it should not be used until it is opened with calling open() function.
		///
		session();
		///
		/// Copy a session object, note - it copies only the reference to the underlying connection, so you should
		/// be very careful when you do it.
		///
		session(session const &);
		///
		/// Assign a session object, note - it copies only the reference to the underlying connection, so you should
		/// be very careful when you do it.
		///
		session const &operator=(session const &);
		///
		/// Destroys the session object, if connection pool is used it returns the object to connection pool.
		///
		/// Note: the connection would not be returned to the pool until all statement and result objects
		/// created using this session are not destroyed.
		///
		~session();

		///
		/// Create a session using a parsed connection string \a ci
		///
		/// \copydetails cppdb::parse_connection_string(std::string const&,std::string&,std::map<std::string,std::string>&);
		///
		session(connection_info const &ci);
		///
		/// Create a session using a connection string \a cs.
		///
		/// \copydetails cppdb::parse_connection_string(std::string const&,std::string&,std::map<std::string,std::string>&);
		///
		session(std::string const &cs);
		///
		/// Create a session using a parsed connection string \a ci and call \a f if
		/// a \ref once() was not called yet.
		///
		/// It is useful for setting session specific options for new
		/// connection, not reused from the pool one.
		///
		/// Requirements: \ref once_functor is an object that can be created from generic
		/// function like object func, such that func(*this) is valid expression
		///
		/// \copydetails cppdb::parse_connection_string(std::string const&,std::string&,std::map<std::string,std::string>&);
		///
		session(connection_info const &ci,once_functor const &f);
		///
		/// Create a session using a connection string \a cs and call \a f if
		/// a \ref once() was not called yet.
		///
		/// It is useful for setting session specific options for new
		/// connection, not reused from the pool one.
		///
		/// Requirements: \ref once_functor is an object that can be created from generic
		/// function like object func, such that func(*this) is valid expression
		///
		/// \copydetails cppdb::parse_connection_string(std::string const&,std::string&,std::map<std::string,std::string>&);
		///
		session(std::string const &cs,once_functor const &f);
		///
		/// Create a session using a pointer to backend::connection and call \a f if
		/// a \ref once() was not called yet.
		///
		/// It is useful for setting session specific options for new
		/// connection, not reused from the pool one.
		///
		/// Requirements: \ref once_functor is an object that can be created from generic
		/// function like object func, such that func(*this) is valid expression
		///
		///
		session(ref_ptr<backend::connection> conn,once_functor const &f);
		///
		/// Create a session using a pointer to backend::connection.
		///
		session(ref_ptr<backend::connection> conn);

		///
		/// Open a session using a connection_info object - parsed connection string \a ci.
		///
		void open(connection_info const &ci);
		///
		/// Open a session using a connection string \a cs.
		///
		/// \copydetails cppdb::parse_connection_string(std::string const&,std::string&,std::map<std::string,std::string>&);
		///
		void open(std::string const &cs);
		///
		/// Close current connection, note, if connection pooling is used the connection is not actually become closed but
		/// rather recycled for future use.
		///
		void close();
		///
		/// Check if the session was opened.
		///
		bool is_open();

		///
		/// Create a new statement, by default is creates prepared statement - create_prepared_statement() unless \@use_prepared connection string
		/// property is set to off, then it uses normal statements by calling create_statement()
		///
		/// This is the most convenient function to create statements with.
		///
		statement prepare(std::string const &query);
		///
		/// Syntactic sugar, same as prepare(q)
		///
		statement operator<<(std::string const &q);
		///
		/// Syntactic sugar, same as prepare(s)
		///
		statement operator<<(char const *s);


		///
		/// Create ordinary statement it generally unprepared statement and it is never cached. It should
		/// be used when such statement is executed rarely or very customized.
		///
		statement create_statement(std::string const &q);
		///
		/// Create prepared statement that will be cached for next calls.
		///
		statement create_prepared_statement(std::string const &q);
		///
		/// Create prepared statement however don't use statements cache and for it. Useful for creation
		/// of custom or rarely executed statements that should be executed several times at this point in program.
		///
		statement create_prepared_uncached_statement(std::string const &q);

		///
		/// Remove all statements from the cache.
		///
		void clear_cache();

		///
		/// Clear connections pool associated with this session's connection.
		///
		/// Automatically calls clear_cache();
		///
		void clear_pool();

		///
		/// Begin a transaction. Don't use it directly for RAII reasons. Use transaction class instead.
		///
		void begin();
		///
		/// Commit a transaction. Don't use it directly for RAII reasons. Use transaction class instead.
		///
		void commit();
		///
		/// Rollback a transaction. Don't use it directly for RAII reasons. Use transaction class instead.
		///
		void rollback();

		///
		/// Escape a string in range [\a b,\a e) for inclusion in SQL statement. It does not add quotation marks at beginning and end.
		/// It is designed to be used with text, don't use it with generic binary data.
		///
		/// Some backends (odbc) may not support this.
		///
		std::string escape(char const *b,char const *e);
		///
		/// Escape a NULL terminated string \a s for inclusion in SQL statement. It does not add quotation marks at beginning and end.
		/// It is designed to be used with text, don't use it with generic binary data.
		///
		/// Some backends (odbc) may not support this.
		///
		std::string escape(char const *s);
		///
		/// Escape a string \a s for inclusion in SQL statement. It does not add quotation marks at beginning and end.
		/// It is designed to be used with text, don't use it with generic binary data.
		///
		/// Some backends (odbc) may not support this.
		///
		std::string escape(std::string const &s);
		///
		/// Get the driver name, as mysql, postgresql, odbc, sqlite3.
		///
		std::string driver();
		///
		/// Get an SQL engine name, it may be not the same as driver name for multiple engine drivers like odbc.
		///
		std::string engine();
		///
		/// Check if this session's connection can be recycled for reuse in a pool.
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

		///
		/// Returns true of session specific initialization is done, otherwise returns false
		///
		bool once_called();
		///
		/// Set flag to true if session specific initialization is done, otherwise set it to false
		///
		void once_called(bool state);
		///
		/// Call the functional \a f on the connection only once. If the connection
		/// created first time then f would be called. If the connection is created
		/// from connection pool and thus setup functor was called, f would not be called.
		///
		/// Requirements: \ref once_functor is an object that can be created from generic
		/// function like object func, such that func(*this) is valid expression
		///
		void once(once_functor const &f);


		///
		/// Get connection specific object by its type \a t, returns 0 if not installed yet
		///
		connection_specific_data *get_specific(std::type_info const &t);
		///
		/// Transfers ownership on the connection specific object of type \a t, returns 0 if not installed yet
		///
		connection_specific_data *release_specific(std::type_info const &t);
		///
		/// Deletes connection specific object of type \a t, and sets a new one \a p (if not NULL)
		///
		void reset_specific(std::type_info const &t,connection_specific_data *p=0);

		///
		/// Get connection specific object by its type \a T, returns 0 if not installed yet
		///
		template<typename T>
		T *get_specific()
		{
			return static_cast<T*>(get_specific(typeid(T)));
		}
		///
		/// Transfers ownership on the connection specific object of type \a T, returns 0 if not installed yet
		///
		template<typename T>
		T *release_specific()
		{
			return static_cast<T*>(release_specific(typeid(T)));
		}
		///
		/// Deletes connection specific object of type \a T, and sets a new one \a p (if not NULL)
		///
		template<typename T>
		void reset_specific(T *p=0)
		{
			reset_specific(typeid(T),p);
		}

	private:
		struct data;
		std::unique_ptr<data> d;
		ref_ptr<backend::connection> conn_;
	};

	///
	/// \brief The transaction guard
	///
	/// This class is RAII transaction guard that causes automatic transaction rollback on stack unwind, unless
	/// the transaction is committed
	///
	class CPPDB_API transaction {
		transaction(transaction const &);
		void operator=(transaction const &);
	public:
		///
		/// Begin a transaction on session \a s, calls s.begin()
		///
		transaction(session &s);
		///
		/// If the transaction wasn't committed or rolled back calls session::rollback() for the session it was created with.
		///
		~transaction();
		///
		/// Commit a transaction on the session.  Calls session::commit() for the session it was created with.
		///
		void commit();
		///
		/// Rollback a transaction on the session.  Calls session::rollback() for the session it was created with.
		///
		void rollback();
	private:

		struct data;
		session *s_;
		bool commited_;
		std::unique_ptr<data> d;
	};


} // cppdb

#endif
