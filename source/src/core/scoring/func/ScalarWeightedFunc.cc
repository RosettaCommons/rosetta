// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/ScalarWeightedFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <sstream>
#include <iostream>
#include <string>


namespace core {
namespace scoring {
namespace func {

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


	void
	ScalarWeightedFunc::show_definition( std::ostream &out ) const
	{
		out << "SCALARWEIGHTEDFUNC  " << weight_ << " ";
		func_to_weight_->show_definition(out);
	}

	Size
	ScalarWeightedFunc::show_violations(std::ostream &out, Real x, Size verbose_level, Real threshold) const
	{
		out << "SCALARWEIGHTEDFUNC with weight:  " << weight_ << std::endl;
		func_to_weight_->show_definition(out);
		out << " with verbose_level " << verbose_level << ", threshold " << threshold << " and weighted_score " << ScalarWeightedFunc::func(x) << std::endl;
		return func_to_weight_->show_violations(out,x,verbose_level,threshold);
	}


} // namespace constraints
} // namespace scoring
} // namespace core
