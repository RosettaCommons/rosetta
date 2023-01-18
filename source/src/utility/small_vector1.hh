// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/small_vector1.hh
/// @brief  boost::container::small_vector with 1-based indexing
/// @author Jack Maguire


#ifndef INCLUDED_utility_small_vector1_HH
#define INCLUDED_utility_small_vector1_HH

#include <boost/container/small_vector.hpp>
#include <utility/assert.hh>

#if defined(__clang__) || defined(__llvm__)
#if __clang_major__ < 4
#include <utility/vector1.hh>
#include <utility/vector0.hh>

//see https://github.com/RosettaCommons/main/pull/5076 for why we do this
#define SMALL_VECTOR_CLANG_HACK

namespace utility {

template< typename T, std::size_t BUFFER_SIZE >
using small_vector0 = utility::vector0< T >;
	
template< typename T, std::size_t BUFFER_SIZE >
using small_vector1 = utility::vector1< T >;
	
} //namespace utility

#endif
#endif

#ifndef SMALL_VECTOR_CLANG_HACK

namespace utility {

template< typename T, std::size_t L, std::size_t BUFFER_SIZE >
class small_vectorL :
    public boost::container::small_vector< T, BUFFER_SIZE >
{
public:
	using boost::container::small_vector< T, BUFFER_SIZE >::small_vector;
	small_vectorL( small_vectorL< T, L, BUFFER_SIZE > const & ) = default;

	using boost::container::small_vector< T, BUFFER_SIZE >::size;

	using boost::container::small_vector< T, BUFFER_SIZE >::operator=;
	small_vectorL & operator=( small_vectorL< T, L, BUFFER_SIZE > const & ) = default;
	    
	T & operator[]( std::size_t const i ){
		debug_assert( i >= L );
		debug_assert( i - L < size() );
		return boost::container::small_vector< T, BUFFER_SIZE >::operator[]( i-L );
	}

	T const & operator[]( std::size_t const i ) const {
		debug_assert( i >= L );
		debug_assert( i - L < size() );
		return boost::container::small_vector< T, BUFFER_SIZE >::operator[]( i-L );
	}
	    
	T & at( std::size_t const i ){
		debug_assert( i >= L );
		debug_assert( i - L < size() );
		return boost::container::small_vector< T, BUFFER_SIZE >::at( i-L );
	}

	T const & at( std::size_t const i ) const {
		debug_assert( i >= L );
		debug_assert( i - L < size() );
		return boost::container::small_vector< T, BUFFER_SIZE >::at( i-L );
	}
};

template< typename T, std::size_t BUFFER_SIZE >
using small_vector0 = small_vectorL< T, 0, BUFFER_SIZE >;
	
template< typename T, std::size_t BUFFER_SIZE >
using small_vector1 = small_vectorL< T, 1, BUFFER_SIZE >;

} // namespace utility

#endif //ndef SMALL_VECTOR_CLANG_HACK

#endif // INCLUDED_utility_vector1_HH
