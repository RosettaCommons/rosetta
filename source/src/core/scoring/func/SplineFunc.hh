// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/func/SplineFunc.hh
/// @brief  Constraint function for looking up data in a histogram over which a spline is created
/// @author Stephanie Hirst (stephanie.j.hirst@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_func_SplineFunc_hh
#define INCLUDED_core_scoring_func_SplineFunc_hh

// Unit Headers
#include <core/scoring/func/SplineFunc.fwd.hh>

// Package Headers
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/Func.hh>

// Project Headers

// Utility and Numeric Headers
#include <utility/vector1.hh>
#include <numeric/interpolation/spline/Interpolator.hh>

// C++ Headers
#include <iostream>


namespace core {
	namespace scoring {
		namespace func {

class SplineFunc : public Func {

public:

	// @brief SplineFunc construction and destruction
	SplineFunc();
	~SplineFunc();

	/// @brief returns a clone of this SplineFunc
	virtual
	FuncOP clone() const { return FuncOP( new SplineFunc( *this ) ); }

	/// @brief return SplineFunc member variables
	core::Real get_exp_val();

	std::string get_filename();

	std::string get_KB_description();

	core::Real get_weight();

	core::Real get_bin_size();

	core::Real get_lower_bound_x();

	core::Real get_upper_bound_x();

	core::Real get_lower_bound_y();

	core::Real get_upper_bound_y();

	core::Real get_lower_bound_dy();

	core::Real get_upper_bound_dy();

	/// @brief initialize this SplineFunc from the given izstream.
	virtual
	void read_data ( std::istream &in );

	/// @brief Returns the value of this SplineFunc evaluated at distance x.
	virtual
	core::Real func( core::Real const x ) const;

	/// @brief Returns the value of the first derivative of this SplineFunc at distance x.
	virtual
	core::Real dfunc( core::Real const x ) const;

	/// @brief show the definition of this SplineFunc to the specified output stream.
	virtual
	void show_definition( std::ostream &out ) const;

	/// @brief show some sort of stringified representation of the violations for this constraint.
	virtual
	core::Size show_violations( std::ostream &out, core::Real x, core::Size verbose_level, core::Real threshold = 1 ) const;

private:

	core::Real exp_val_;
	std::string filename_;
	std::string KB_description_;
	core::Real weight_;
	core::Real bin_size_;
	core::Real lower_bound_x_;
	core::Real lower_bound_y_;
	core::Real upper_bound_x_;
	core::Real upper_bound_y_;
	core::Real lower_bound_dy_;
	core::Real upper_bound_dy_;
	utility::vector1<core::Real> bins_vect_;
	utility::vector1<core::Real>::size_type bins_vect_size_;
	utility::vector1<core::Real> potential_vect_;
	utility::vector1<core::Real>::size_type potential_vect_size_;
	numeric::interpolation::spline::InterpolatorOP interpolator_;

}; // SplineFunc class

		} // constraints
	} // scoring
} // core

#endif
