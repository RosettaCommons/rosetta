// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularGeneral1D_Func.hh
/// @brief A general 1D function that can be initialized by FArray or from text file. There's stuff like this all over Rosetta.
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_func_CircularGeneral1D_Func_HH
#define INCLUDED_core_scoring_func_CircularGeneral1D_Func_HH

#include <core/scoring/func/Func.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

/// @brief Function that allows return of arbitrary FArrays -- this time circularized.
class CircularGeneral1D_Func : public Func {
public:
	CircularGeneral1D_Func(
		ObjexxFCL::FArray1D< core::Real > const & data,
		core::Real const & xmin,
		core::Real const & xbin
	);

	CircularGeneral1D_Func( std::string const & filename );

	FuncOP clone() const { return FuncOP( new CircularGeneral1D_Func( *this ) ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

private:

	ObjexxFCL::FArray1D< core::Real > data_;
	Real xmin_, xbin_;
	Size num_bins_;

};

} // constraints
} // scoring
} // core

#endif
