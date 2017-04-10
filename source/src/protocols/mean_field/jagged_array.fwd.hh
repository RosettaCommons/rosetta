// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/mean_field/jagged_array.fwd.hh
/// @brief  protocols::mean_field::jagged_array forward declarations
/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)


#ifndef INCLUDED_mean_field_jagged_array_FWD_HH
#define INCLUDED_mean_field_jagged_array_FWD_HH


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


namespace protocols {
namespace mean_field {


// Forward
template< typename T, typename A = std::allocator< T > > class jagged_array;


// Types
typedef  jagged_array < bool >               jagged_array_bool;
typedef  jagged_array< short int >           jagged_array_short;
typedef  jagged_array< int >                 jagged_array_int;
typedef  jagged_array< long int >            jagged_array_long;
typedef  jagged_array< unsigned short int >  jagged_array_ushort;
typedef  jagged_array< unsigned int >        jagged_array_uint;
typedef  jagged_array< unsigned long int >   jagged_array_ulong;
typedef  jagged_array< std::size_t >         jagged_array_size_t;
typedef  jagged_array< std::size_t >         jagged_array_size;
typedef  jagged_array< ssize_t >             jagged_array_ssize_t;
typedef  jagged_array< ssize_t >             jagged_array_ssize;
typedef  jagged_array< float >               jagged_array_float;
typedef  jagged_array< double >              jagged_array_double;
typedef  jagged_array< long double >         jagged_array_longdouble;
typedef  jagged_array< char >                jagged_array_char;
typedef  jagged_array< unsigned char >       jagged_array_uchar;
typedef  jagged_array< signed char >         jagged_array_schar;

} // namespace mean_field
} // namespace protocols

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#endif

#endif // INCLUDED_mean_field_jagged_array_FWD_HH
