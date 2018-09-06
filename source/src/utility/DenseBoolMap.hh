// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/DenseBoolMap.hh
/// @brief Provides a memory efficient way to store a large number of boolean values.
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_utility_dense_bool_map_hh
#define INCLUDED_utility_dense_bool_map_hh

#include <array>
#include <utility/exit.hh>

namespace utility {

///@brief This struct condenses N boolean values into roughly N bits, saving roughly 4x memory.
/// NUM_ELEMS parameter takes in the number of boolean values.
/// BASE_INDEX parameter takes in the index of the first value.
/// The recommended way to use this struct is with an enum. See core::scoring::hbonds::graph::AtomInfo for an example
template < unsigned int NUM_ELEMS, unsigned int BASE_INDEX >
struct DenseBoolMap
{
public:
	DenseBoolMap(){
#ifndef NDEBUG
		if ( num_bytes() != data_.size() ) {
			utility_exit_with_message( "Discrepancy between num_bytes() and data_.size()" );
		}
#endif
		for ( unsigned int i = 0; i < num_bytes(); ++i ) {
			data_[ i ] = 0;
		}
	}

	static constexpr unsigned int num_bytes(){
		return ( NUM_ELEMS + 3 ) / 4;
	}

	static constexpr unsigned int byte_for_element( unsigned int const element ){
		return ( element - BASE_INDEX ) / 4;
	}

	static constexpr unsigned char mask_for_element( unsigned int const element ){
		return 1 << ( ( element - BASE_INDEX ) % 4 );
	}

	inline void set( unsigned int const element, bool const setting ){
#ifndef NDEBUG
		if ( element - BASE_INDEX >= NUM_ELEMS ) {
			utility_exit_with_message( "DenseBoolMap element " + std::to_string( element ) + " is out of range." );
		}
#endif
		if ( setting ) {
			data_[ byte_for_element( element ) ] |= mask_for_element( element );
		} else {
			data_[ byte_for_element( element ) ] &= ( mask_for_element( element ) ^ 0xF );
		}
	}

	///@brief templated equivalent to set( element, setting ). I personally like this one better
	/// because it does bounds-checking at compile time and not just in debug mode. Of course,
	/// this only works if you know the value for "element" when calling set().
	template< unsigned int const element >
	inline void set( bool const setting ){
		static_assert( element - BASE_INDEX < NUM_ELEMS, "element is not in range for DenseBoolMap" );
		if ( setting ) {
			data_[ byte_for_element( element ) ] |= mask_for_element( element );
		} else {
			data_[ byte_for_element( element ) ] &= ( mask_for_element( element ) ^ 0xF );
		}
	}

	inline bool get( unsigned int const element ) const {
#ifndef NDEBUG
		if ( element - BASE_INDEX >= NUM_ELEMS ) {
			utility_exit_with_message( "DenseBoolMap element " + std::to_string( element ) + " is out of range." );
		}
#endif
		return data_[ byte_for_element( element ) ] & mask_for_element( element );
	}

	///@brief templated equivalent to get( element ). I personally like this one better because
	/// it does bounds-checking at compile time and not just in debug mode. Of course, this only
	/// works if you know the value for "element" when calling get().
	template< unsigned int const element >
	constexpr bool get() const {
		static_assert( element - BASE_INDEX < NUM_ELEMS, "element is not in range for DenseBoolMap" );
		return data_[ byte_for_element( element ) ] & mask_for_element( element );
	}


private:
	//Please change this line back to the following once all compilers can handle constexpr:
	//std::array< unsigned char, num_bytes() > data_;
	std::array< unsigned char, ( NUM_ELEMS + 3 ) / 4 > data_;

};

}

#endif
