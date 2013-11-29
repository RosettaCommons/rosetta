// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SkipViolFunc.cc
/// @brief Weighted constraint function that reweights other constraints
/// by a constant scalar value.
/// @author James Thompson, Greg Taylor


#include <core/scoring/func/SkipViolFunc.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers
#include <string>

namespace core {
namespace scoring {
namespace constraints {

	Real
	SkipViolFunc::func( Real const x ) const
	{
		return func_to_weight_->func( x ) * weight_;
	}

	Real
	SkipViolFunc::dfunc( Real const x ) const
	{
		return func_to_weight_->dfunc( x ) * weight_;
	}

	void
	SkipViolFunc::read_data( std::istream& in )
	{
		in >> count_viols_;
		weight_ = 1.0;
		FuncFactory func_factory;
		std::string func_type;
		in >> func_type;
		func_to_weight_ = func_factory.func_types_[ func_type ]->clone();
    func_to_weight_->read_data( in );
	}

   /// @brief show some sort of stringified representation of the violations for this constraint.
core::Size SkipViolFunc::show_violations(
			std::ostream& out,
			Real r,
			Size verbose_level,
			Real threshold
) const {
	Size ct ( func_to_weight_->show_violations( out, r, verbose_level, threshold ) );
	return ct;
}

void
SkipViolFunc::show_definition( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	out << "SKIPVIOLFUNC" << RJ( 7, count_viols_ ) << " ";
	func_to_weight_->show_definition( out );
}

} // namespace constraints
} // namespace scoring
} // namespace core

