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
#include <cppdb/shared_object.h>
#include <string>
#include <utility>

#if defined(CPPDB_DISABLE_SHARED_OBJECT_LOADING)
// Disable shared object loading by simply asserting false when
// an attempt to call dlopen, dlclose or dlsym is made
#	define RTLD_LAZY 0
#	include <cassert>
	namespace cppdb {
		namespace {
                	void *dlopen(char const *,int)
			{
                        	assert(false);
                        	char * return_value = NULL;
				return return_value; // have to return something.
			}
                	void dlclose(void *)
			{
                        	assert(false);
			}
                	void *dlsym(void *,char const *)
			{
                        	assert(false); 
                        	char * return_value = NULL;
				return return_value; // have to return something.
			}
		}
	}

#elif defined(WIN32) || defined(_WIN32) || defined(__WIN32)
#	include <windows.h>
#	define RTLD_LAZY 0

	namespace cppdb {
		namespace {
			void *dlopen(char const *name,int /*unused*/)
			{
				return LoadLibrary(name);
			}
			void dlclose(void *h)
			{
				HMODULE m=(HMODULE)(h);
				FreeLibrary(m);
			}
			void *dlsym(void *h,char const *sym)
			{
				HMODULE m=(HMODULE)(h);
				return (void *)GetProcAddress(m,sym);
			}
		}
	}

#else
#	include <dlfcn.h>
#endif
namespace cppdb {
	shared_object::shared_object(std::string const & name,void *h) :
		dlname_(name),
		handle_(h)
	{
	}
	shared_object::~shared_object()
	{
		dlclose(handle_);
	}
	ref_ptr<shared_object> shared_object::open(std::string const &name)
	{
		ref_ptr<shared_object> dl;
		void *h=dlopen(name.c_str(),RTLD_LAZY);
		if(!h) {
			return dl;
		}
		try {
			dl.reset(new shared_object(name,h));
			h=nullptr;
			return dl;
		}
		catch(...) {
			if(h) {
				dlclose(h);
			}
			throw;
		}
	}
	void *shared_object::sym(std::string const &s)
	{
		return dlsym(handle_,s.c_str());
	}
	void *shared_object::safe_sym(std::string const &s)
	{
		void *p=sym(s);
		if(!p) {
			throw cppdb_error("cppdb::shared_object::failed to resolve symbol [" + s +"] in " + dlname_);
		}
		return p;
	}
} // cppdb
