// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/internal/RowsPointer.hh
/// @brief  Contiguous row-ordered 3x3 matrix pointer wrapper class
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_internal_RowsPointer_hh
#define INCLUDED_numeric_internal_RowsPointer_hh


// Package headers
#include <numeric/xyzMatrix.fwd.hh>


namespace numeric {


template< typename T >
class RowsPointer
{


private: // Friends


	template< typename > friend class numeric::xyzMatrix;


public: // Creation


	/// @brief Contiguous row-ordered 3x3 xyzMatrix pointer constructor
	/// @warning No way to check that argument points to nine values
	inline
	explicit
	RowsPointer( T const * p_a ) :
		p_( p_a )
	{}


public: // Methods: conversion operators


	/// @brief Conversion to wrapped pointer
	inline
	operator T const *() const
	{
		return p_;
	}


private: // Fields


	/// @brief Pointer (non-owning) to contiguous row-ordered 3x3 matrix
	T const * p_;


}; // RowsPointer


} // namespace numeric


#endif // INCLUDED_numeric_internal_RowsPointer_HH
