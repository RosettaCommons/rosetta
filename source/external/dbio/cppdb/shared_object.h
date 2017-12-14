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
#ifndef CPPDB_SHARED_OBJECT_H
#define CPPDB_SHARED_OBJECT_H

#include <cppdb/defs.h>
#include <cppdb/ref_ptr.h>


namespace cppdb {
	///
	/// \brief This class allows to load and unload shared objects in simple and exception safe way.
	///
	class CPPDB_API shared_object : public ref_counted {
		shared_object() : handle_(0) {}
		shared_object(std::string const & name,void *h);
		shared_object(shared_object const &);
		void operator=(shared_object const &);
	public:
		~shared_object();
		///
		/// Load shared object, returns empty pointer if the object does not exits or not loadable
		///
		static ref_ptr<shared_object> open(std::string const &name);
		///
		/// Resolve symbol \a name and return pointer on it, throws cppdb_error if the symbol can't be resolved 
		///
		void *safe_sym(std::string const &name);

		///
		/// Resolve symbol \a name and return pointer on it, returns NULL if the symbol can't be resolved 
		///
		void *sym(std::string const &name);

		///
		/// Resolve symbol \a name and assign it to \a v, returns false if the symbol can't be resolved 
		///
		template<typename T>
		bool resolve(std::string const &s,T *&v)
		{
			void *p=sym(s);
			if(!p) {
				return false;
			}
			v=(T*)(p);
			return true;
		}
		///
		/// Resolve symbol \a name and assign it to v, throws cppdb_error if the symbol can't be resolved 
		///
		template<typename T>
		void safe_resolve(std::string const &s,T *&v)
		{
			v=(T*)(sym(s));
		}

	private:
		std::string dlname_;
		void *handle_;
	};
}


#endif
