// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/constraints/ConstantConstraint.cc

#include <core/types.hh>
#include <core/scoring/constraints/ConstantConstraint.hh>
#include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <core/id/AtomID.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

/// @brief compute score
Real
ConstantConstraint::score() const {
	return func_->func(0);
}

std::string ConstantConstraint::type() const {
	return "Constant";
}

/// @brief compute score
void
ConstantConstraint::score( XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += score();
}

/// @brief compute atom deriv
void
ConstantConstraint::fill_f1_f2(
	AtomID const & ,
	XYZ_Func const &,
	Vector & /*F1*/,
	Vector & /*F2*/,
	EnergyMap const &
) const {}

/// @brief Constructor
ConstantConstraint::ConstantConstraint(
	FuncOP func_in, // we take ownership of this guy
	ScoreType scotype
):
	Constraint( scotype ),
	func_( func_in )
{}

ConstantConstraint::~ConstantConstraint() {}

/// @brief number of atoms --- zero
Size
ConstantConstraint::natoms() const
{
	return 0;
}

id::AtomID const &
ConstantConstraint::atom( core::Size const /*n*/ ) const {
	utility_exit_with_message( "ConstantConstraint::atom() - no atoms exist" );
	return *(new AtomID); // never happens
}

/// @brief output violation of constraint (none!)
Size ConstantConstraint::show_violations( std::ostream &, pose::Pose const &, Size, Real  ) const {
	return 0;
}

void ConstantConstraint::show( std::ostream& out ) const {
	out << "ConstantConstraint" << std::endl;
	func_->show( out );
}

} // constraints
} // scoring
} // core
