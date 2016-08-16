// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/internal/RowPointers.hh
/// @brief  3x3 matrix row pointers wrapper class
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_internal_RowPointers_hh
#define INCLUDED_numeric_internal_RowPointers_hh


// Package headers
#include <numeric/xyzMatrix.fwd.hh>


namespace numeric {


template< typename T >
class RowPointers
{


private: // Friends


	template< typename > friend class numeric::xyzMatrix;


public: // Creation


	/// @brief Row pointers constructor
	/// @warning No way to check that arguments each point to three values
	inline
	RowPointers(
		T const * xp_a, // Pointer to x row
		T const * yp_a, // Pointer to y row
		T const * zp_a  // Pointer to z row
	) :
		xp_( xp_a ),
		yp_( yp_a ),
		zp_( zp_a )
	{}


public: // Properties


	/// @brief x row pointer
	inline
	T const *
	xp() const
	{
		return xp_;
	}


	/// @brief y row pointer
	inline
	T const *
	yp() const
	{
		return yp_;
	}


	/// @brief z row pointer
	inline
	T const *
	zp() const
	{
		return zp_;
	}


private: // Fields


	/// @brief Pointer (non-owning) to x row
	T const * xp_;

	/// @brief Pointer (non-owning) to y row
	T const * yp_;

	/// @brief Pointer (non-owning) to z row
	T const * zp_;


}; // RowPointers


} // namespace numeric


#endif // INCLUDED_numeric_internal_RowPointers_HH
