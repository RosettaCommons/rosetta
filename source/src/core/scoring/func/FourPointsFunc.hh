// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_func_FourPointsFunc_hh
#define INCLUDED_core_scoring_func_FourPointsFunc_hh

// Unit headers
#include <core/scoring/func/FourPointsFunc.fwd.hh>

// Package headers
#include <core/scoring/func/XYZ_Func.hh>

// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>

// Numeric Headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>

#include <utility/vector1_bool.hh>



namespace core {
namespace scoring {
namespace func {

/// @brief A simple class that represents the coordinates of four points, pretending
/// that they all belong to residue 1.  The residue() method is not implemented and
/// cause a utility_exit.
class FourPointsFunc : public func::XYZ_Func {
public:
	typedef func::XYZ_Func parent;
	typedef parent::AtomID AtomID;
	typedef parent::Residue Residue;
	typedef parent::Conformation Conformation;

public:
	FourPointsFunc();

	virtual
	~FourPointsFunc();

	/// @brief set the coordinate for one of the four atoms
	void xyz( Size atomid, Vector const & coord );

	virtual
	Vector const &
	operator()( AtomID const & id ) const;


	virtual
	Residue const &
	residue( Size seqpos ) const;

private:
	utility::vector1< Vector > points_;
};



} // constraints
} // scoring
} // core

#endif
