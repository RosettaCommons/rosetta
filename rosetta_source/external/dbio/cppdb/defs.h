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
#ifndef CPPDB_DEFS_H
#define CPPDB_DEFS_H

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(__CYGWIN__)
#  if defined(DLL_EXPORT) || defined(CPPDB_EXPORTS) || defined(CPPDB_DRIVER_EXPORTS)
#    ifdef CPPDB_SOURCE
#      define CPPDB_API __declspec(dllexport)
#    else
#      define CPPDB_API __declspec(dllimport)
#    endif
#  endif
#  if defined(DLL_EXPORT) || defined(CPPDB_DRIVER_EXPORTS)
#    ifdef CPPDB_DRIVER_SOURCE
#      define CPPDB_DRIVER_API __declspec(dllexport)
#    else
#      define CPPDB_DRIVER_API __declspec(dllimport)
#    endif
#  endif
#endif


#ifndef CPPDB_API
#  define CPPDB_API
#endif

#ifndef CPPDB_DRIVER_API
#  define CPPDB_DRIVER_API
#endif

#endif
