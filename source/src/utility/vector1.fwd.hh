// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vector1.fwd.hh
/// @brief  utility::vector1 forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vector1_FWD_HH
#define INCLUDED_utility_vector1_FWD_HH


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
template< typename T, typename A = std::allocator< T > > class vector1;


// Types
typedef  vector1< bool >                vector1_bool;
typedef  vector1< short int >           vector1_short;
typedef  vector1< int >                 vector1_int;
typedef  vector1< long int >            vector1_long;
typedef  vector1< unsigned short int >  vector1_ushort;
typedef  vector1< unsigned int >        vector1_uint;
typedef  vector1< unsigned long int >   vector1_ulong;
typedef  vector1< std::size_t >         vector1_size_t;
typedef  vector1< std::size_t >         vector1_size;
typedef  vector1< platform::SSize >     vector1_ssize_t;
typedef  vector1< platform::SSize >     vector1_ssize;
typedef  vector1< float >               vector1_float;
typedef  vector1< double >              vector1_double;
typedef  vector1< long double >         vector1_longdouble;
typedef  vector1< char >                vector1_char;
typedef  vector1< unsigned char >       vector1_uchar;
typedef  vector1< signed char >         vector1_schar;


} // namespace utility

#endif // INCLUDED_utility_vector1_FWD_HH
