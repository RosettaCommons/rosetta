// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/conversions.hh
/// @brief  Conversions between degrees and radians
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_conversions_hh
#define INCLUDED_numeric_conversions_hh


// Package headers
#include <numeric/NumericTraits.hh>


namespace numeric {
namespace conversions {

enum AngleUnit { RADIANS, DEGREES };

/// @brief Radians of degrees
template< typename T >
inline
T
radians( T const & degrees )
{
	return degrees * NumericTraits< T >::degrees_to_radians();
}


/// @brief Radians of any angle
template< typename T >
inline
T
radians( T const & angle, AngleUnit const unit )
{
	if ( unit == DEGREES ) {
		return angle * NumericTraits< T >::degrees_to_radians();
	} else /* unit == RADIANS */ {
		return angle;
	}
}


/// @brief Radians from degrees
template< typename T >
inline
T &
to_radians( T & degrees )
{
	return degrees *= NumericTraits< T >::degrees_to_radians();
}


/// @brief Radians from any angle
template< typename T >
inline
T &
to_radians( T & angle, AngleUnit const unit )
{
	if ( unit == DEGREES ) {
		return angle *= NumericTraits< T >::degrees_to_radians();
	} else /* unit == RADIANS */ {
		return angle;
	}
}


/// @brief Degrees of radians
template< typename T >
inline
T
degrees( T const & radians )
{
	return radians * NumericTraits< T >::radians_to_degrees();
}


/// @brief Degrees of any angle
template< typename T >
inline
T
degrees( T const & angle, AngleUnit const unit )
{
	if ( unit == RADIANS ) {
		return angle * NumericTraits< T >::radians_to_degrees();
	} else /* unit == DEGREES */ {
		return angle;
	}
}


/// @brief Degrees from radians
template< typename T >
inline
T &
to_degrees( T & radians )
{
	return radians *= NumericTraits< T >::radians_to_degrees();
}


/// @brief Degrees from any angle
template< typename T >
inline
T &
to_degrees( T & angle, AngleUnit const unit )
{
	if ( unit == RADIANS ) {
		return angle *= NumericTraits< T >::radians_to_degrees();
	} else /* unit == DEGREES */ {
		return angle;
	}
}


/// @brief Any angle from radians
template< typename T >
inline
T
from_radians( T const & angle, AngleUnit const unit )
{
	if ( unit == DEGREES ) {
		return angle * NumericTraits< T >::radians_to_degrees();
	} else /* unit == RADIANS */ {
		return angle;
	}
}


/// @brief Any angle from radians
template< typename T >
inline
T
from_degrees( T const & angle, AngleUnit const unit )
{
	if ( unit == RADIANS ) {
		return angle * NumericTraits< T >::degrees_to_radians();
	} else /* unit == DEGREES */ {
		return angle;
	}
}


} // namespace conversions
} // namespace numeric


#endif // INCLUDED_numeric_conversions_HH
