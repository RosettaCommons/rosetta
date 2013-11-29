// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/ConstantFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author John Karanicolas

#include <core/scoring/func/ConstantFunc.hh>
#include <core/types.hh>

// C++ Headers
#include <iostream>

namespace core {
namespace scoring {
namespace constraints {

	void
	ConstantFunc::read_data( std::istream & in ) {
		in 	>> return_val_;
	}

	Real
	ConstantFunc::func( Real const ) const	{
		return return_val_;
	} // func

	Real
	ConstantFunc::dfunc( Real const ) const {
		return 0;
	} // dfunc

	void ConstantFunc::show_definition( std::ostream & out ) const {
		out << "CONSTANTFUNC " << return_val_ << "\n";
	}

} // namespace constraints
} // namespace scoring
} // namespace core

