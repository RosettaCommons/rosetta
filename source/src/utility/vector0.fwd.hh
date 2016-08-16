// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector0.fwd.hh
/// @brief  utility::vector0 forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vector0_FWD_HH
#define INCLUDED_utility_vector0_FWD_HH


// Platform headers
#include <platform/types.hh>

// C++ headers
#include <cstddef>


// std::allocator Declaration
#ifdef UNUSUAL_ALLOCATOR_DECLARATION

// C++ headers
#include <vector>

#else

// Faster but not 100% portable
namespace std { template< typename > class allocator; }

#endif // UNUSUAL_ALLOCATOR_DECLARATION


namespace utility {


// Forward
template< typename T, typename A = std::allocator< T > > class vector0;


// Types
typedef  vector0< bool >                vector0_bool;
typedef  vector0< short int >           vector0_short;
typedef  vector0< int >                 vector0_int;
typedef  vector0< long int >            vector0_long;
typedef  vector0< unsigned short int >  vector0_ushort;
typedef  vector0< unsigned int >        vector0_uint;
typedef  vector0< unsigned long int >   vector0_ulong;
typedef  vector0< platform::Size >      vector0_size_t;
typedef  vector0< platform::Size >      vector0_size;
typedef  vector0< platform::SSize >     vector0_ssize_t;
typedef  vector0< platform::SSize >     vector0_ssize;
typedef  vector0< float >               vector0_float;
typedef  vector0< double >              vector0_double;
typedef  vector0< long double >         vector0_longdouble;
typedef  vector0< char >                vector0_char;
typedef  vector0< unsigned char >       vector0_uchar;
typedef  vector0< signed char >         vector0_schar;


} // namespace utility


#endif // INCLUDED_utility_vector0_FWD_HH
