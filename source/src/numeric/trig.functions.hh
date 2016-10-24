// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/trig.functions.hh
/// @brief  Trigonometric functions
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_trig_functions_hh
#define INCLUDED_numeric_trig_functions_hh


// Package headers
#include <numeric/numeric.functions.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

//Exception for NaN in sin/cos
#include <utility/excn/Exceptions.hh>


namespace numeric {


/// @brief Secant
template< typename T >
inline
T
sec( T const & x )
{
	return ( T( 1.0 ) / std::cos( x ) );
}


/// @brief Cosecant
template< typename T >
inline
T
csc( T const & x )
{
	return ( T( 1.0 ) / std::sin( x ) );
}


/// @brief Cotangent
template< typename T >
inline
T
cot( T const & x )
{
	return ( T( 1.0 ) / std::tan( x ) );
}


/// @brief Is a sine or cosine value within a specified tolerance of the valid [-1,1] range?
template< typename T >
inline
bool
in_sin_cos_range( T const & x, T const & tol = T( .001 ) )
{
	assert( tol >= T( 0.0 ) );
	return ( ( x >= -( T( 1.0 ) + tol ) ) && ( x <= T( 1.0 ) + tol ) );
}


/// @brief Adjust a sine or cosine value to the valid [-1,1] range if within a specified
///        tolerance or exit with an error.
///
/// @note  The = on the first <= and >= else if conditions were added to work-around
///        optimization register passing bugs where x can be passed by a register
///        with a float value of +/-1.0f but where x has a full register precision
///        value slightly different from +/-1.0f/  The first else if cases were
///        failing but with the = added these if conditions will succeed.
///        The alternative fix of declaring "volatile float const x" is slower for
///        most calls.
///        DON'T REMOVE THESE = EVEN THOUGH THEY APPEAR TO BE SUPERFLUOUS!!!
template< typename T >
inline
T
sin_cos_range( T const & x, T const & tol = T( .001 ) )
{
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::setprecision;
	using std::showpoint;

	assert( tol >= T( 0.0 ) );
	if ( ( x >= T( -1.0 ) ) && ( x <= T( 1.0 ) ) ) { // In valid [-1,+1] range
		return x;
	} else if ( ( x <= T( -1.0 ) ) && ( x >= -( T( 1.0 ) + tol ) ) ) { // Within tolerance
		return T( -1.0 ); // Adjusted value
	} else if ( ( x >= T( 1.0 ) ) && ( x <= T( 1.0 ) + tol ) ) { // Within tolerance
		return T( 1.0 ); // Adjusted value
	} else { // Out of range
		cout << "sin_cos_range ERROR: " << setprecision( 8 ) << showpoint << x << " is outside of [-1,+1] sin and cos value legal range" << endl;
		cerr << "sin_cos_range ERROR: " << setprecision( 8 ) << showpoint << x << " is outside of [-1,+1] sin and cos value legal range" << endl;
#ifdef BOINC
		// There are enough sufficiently weird machines on BOINC that
		// the error was getting triggered fairly often (~5% of jumping runs).
		//SGM This deserves further investigation: Either there is a bug
		//    or we need a larger tolerance for some call sites on certain (which?) h/w
		return ( x >= T( 0.0 ) ? T( 1.0 ) : T( -1.0 ) );
#endif
#ifdef USEMPI
		std::string const warning( "NANs occured in sin_cos_range!" );
		throw( utility::excn::EXCN_Msg_Exception( warning ) );
#endif
		utility_exit();
		return T( 0.0 ); // Keep compiler happy
	}
}

/// like std::acos but with range checking
template< typename T >
inline
T
arccos( T const x )
{
	return std::acos( sin_cos_range( x ) );
}


} // namespace numeric


#endif // INCLUDED_numeric_trig_functions_HH
