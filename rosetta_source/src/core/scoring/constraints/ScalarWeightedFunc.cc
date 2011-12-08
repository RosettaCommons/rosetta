// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/ScalarWeightedFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/constraints/ScalarWeightedFunc.hh>
#include <core/scoring/constraints/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <sstream>


namespace core {
namespace scoring {
namespace constraints {

	ScalarWeightedFunc::~ScalarWeightedFunc() {}

	Real
	ScalarWeightedFunc::func( Real const x ) const
	{
		return func_to_weight_->func( x ) * weight_;
	}



	Real
	ScalarWeightedFunc::dfunc( Real const x ) const
	{
		return func_to_weight_->dfunc( x ) * weight_;
	}

	void
	ScalarWeightedFunc::read_data( std::istream& in )
	{
		in >> weight_;

		FuncFactory func_factory;
		std::string func_type;
		in >> func_type;
		func_to_weight_ = func_factory.func_types_[ func_type ]->clone();
		func_to_weight_->read_data( in );
	}

} // namespace constraints
} // namespace scoring
} // namespace core

