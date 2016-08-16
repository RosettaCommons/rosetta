// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_constraints_DOF_Constraint_hh
#define INCLUDED_core_scoring_constraints_DOF_Constraint_hh

// Unit header
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/DOF_ID.hh>

//Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


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
/// DOF_Constraints are currently unsupported -- they are never evaluated if you
/// put them into a Pose.
class DOF_Constraint : public utility::pointer::ReferenceCount {

public:
	DOF_Constraint(
		id::DOF_ID const & id,
		core::scoring::func::FuncOP func,
		ScoreType t = dof_constraint );

	virtual
	~DOF_Constraint();

	/// @brief Returns the ScoreType
	id::DOF_ID const &
	dof_id() const;

	/// @brief Returns the ScoreType
	ScoreType const &
	score_type() const;

	/// @brief Returns the func::Function
	Real
	func( Real const val ) const;

	/// @brief Returns the func::Function Derivative
	Real
	dfunc( Real const val ) const;

	virtual void show( std::ostream& out ) const;

private:

	id::DOF_ID const dof_id_;
	func::FuncOP func_;
	ScoreType const score_type_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DOF_Constraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_DOF_Constraint )
#endif // SERIALIZATION

#endif
