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


#ifndef INCLUDED_core_scoring_constraints_DOF_Constraint_hh
#define INCLUDED_core_scoring_constraints_DOF_Constraint_hh

// Unit header
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/DOF_ID.hh>

//Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace constraints {



/// @brief This isn't quite a standard constraint since it acts on DOF's directly
///          rather than on XYZ coordinates.
/// @details All DOF_Constraints are expected to be immutable once created,
/// meaning their internal data (state) should not change over their lifetime.
/// This allows DOF_Constraints to be shared between copies of Poses (e.g. in Monte Carlo),
/// and is important for both speed (with thousands of contraints) and correctness.
///
/// To "change" a constraint, remove the old one and add a new and different one.
/// The clone() and steal() methods have been removed because they are
/// unneccessary and incompatible (respectively) with the idea of immutable constraints.
class DOF_Constraint : public utility::pointer::ReferenceCount {

public:
	DOF_Constraint( id::DOF_ID const & id, FuncOP func, ScoreType t = dof_constraint ):
		dof_id_( id ),
		func_( func),
		score_type_(t) {}

	/// @brief Returns the ScoreType
	id::DOF_ID const &
	dof_id() const
	{
		return dof_id_;
	}

	/// @brief Returns the ScoreType
	ScoreType const &
	score_type() const
	{
		return score_type_;
	}

	/// @brief Returns the Function
	Real
	func( Real const val ) const
	{
		return func_->func( val );
	}

	/// @brief Returns the Function Derivative
	Real
	dfunc( Real const val ) const
	{
		return func_->dfunc( val );
	}


	virtual
	~DOF_Constraint(){}

	virtual void show( std::ostream& out ) const {
		out << "DOF_Constraint::show stubbed out!" << std::endl;
	}

private:

	id::DOF_ID const dof_id_;
	constraints::FuncOP func_;
	ScoreType const score_type_;

};


} // constraints
} // scoring
} // core

#endif
