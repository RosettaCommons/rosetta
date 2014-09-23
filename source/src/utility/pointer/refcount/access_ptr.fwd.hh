// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/access_ptr.fwd.hh
/// @brief  utility::pointer::access_ptr forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_pointer_refcount_access_ptr_fwd_hh
#define INCLUDED_utility_pointer_refcount_access_ptr_fwd_hh


// C++ headers
#include <cstddef>


namespace utility {
namespace pointer {


// Forward
template< typename T > class access_ptr;


// Types
typedef  access_ptr< bool >                access_ptr_bool;
typedef  access_ptr< short int >           access_ptr_short;
typedef  access_ptr< int >                 access_ptr_int;
typedef  access_ptr< long int >            access_ptr_long;
typedef  access_ptr< unsigned short int >  access_ptr_ushort;
typedef  access_ptr< unsigned int >        access_ptr_uint;
typedef  access_ptr< unsigned long int >   access_ptr_ulong;
typedef  access_ptr< std::size_t >         access_ptr_size_t;
typedef  access_ptr< float >               access_ptr_float;
typedef  access_ptr< double >              access_ptr_double;
typedef  access_ptr< long double >         access_ptr_longdouble;
typedef  access_ptr< char >                access_ptr_char;
typedef  access_ptr< unsigned char >       access_ptr_uchar;
typedef  access_ptr< signed char >         access_ptr_schar;


} // namespace pointer
} // namespace utility


#endif // INCLUDED_utility_pointer_refcount_access_ptr_FWD_HH
