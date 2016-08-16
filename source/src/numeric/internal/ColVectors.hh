// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/internal/ColVectors.hh
/// @brief  3x3 matrix column vectors wrapper class
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_internal_ColVectors_hh
#define INCLUDED_numeric_internal_ColVectors_hh


// Package headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.hh>


namespace numeric {


template< typename T >
class ColVectors
{


private: // Friends


	template< typename > friend class numeric::xyzMatrix;


public: // Creation


	/// @brief Column vectors constructor
	inline
	ColVectors(
		xyzVector< T > const & x_a, // x column
		xyzVector< T > const & y_a, // y column
		xyzVector< T > const & z_a  // z column
	) :
		x_( x_a ),
		y_( y_a ),
		z_( z_a )
	{}


public: // Properties


	/// @brief x column
	inline
	xyzVector< T > const
	x() const
	{
		return x_;
	}


	/// @brief y column
	inline
	xyzVector< T > const
	y() const
	{
		return y_;
	}


	/// @brief z column
	inline
	xyzVector< T > const
	z() const
	{
		return z_;
	}


private: // Fields


	/// @brief x column
	xyzVector< T > const & x_;

	/// @brief y column
	xyzVector< T > const & y_;

	/// @brief z column
	xyzVector< T > const & z_;


}; // ColVectors


} // namespace numeric


#endif // INCLUDED_numeric_internal_ColVectors_HH
