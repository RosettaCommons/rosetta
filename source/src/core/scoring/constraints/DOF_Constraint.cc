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

// Unit headers
#include <core/scoring/constraints/DOF_Constraint.hh>

// Package headers
#include <core/scoring/func/Func.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

DOF_Constraint::DOF_Constraint(
	id::DOF_ID const & id,
	core::scoring::func::FuncOP func,
	ScoreType t
) :
	dof_id_( id ),
	func_( func),
	score_type_(t)
{}

DOF_Constraint::~DOF_Constraint(){}


/// @brief Returns the ScoreType
id::DOF_ID const &
DOF_Constraint::dof_id() const
{
	return dof_id_;
}

ScoreType const &
DOF_Constraint::score_type() const
{
	return score_type_;
}

/// @brief Returns the func::Function
Real
DOF_Constraint::func( Real const val ) const
{
	return func_->func( val );
}

/// @brief Returns the func::Function Derivative
Real
DOF_Constraint::dfunc( Real const val ) const
{
	return func_->dfunc( val );
}


void DOF_Constraint::show( std::ostream& out ) const {
	out << "DOF_Constraint::show stubbed out!" << std::endl;
}



}
}
}

#ifdef    SERIALIZATION

/// @brief Default constructor that should be used only by the cereal library to deserialize a DOF constraint.
core::scoring::constraints::DOF_Constraint::DOF_Constraint() : dof_id_(), score_type_( dof_constraint ) {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::DOF_Constraint::save( Archive & arc ) const {
	arc( CEREAL_NVP( dof_id_ ) ); // const id::DOF_ID
	arc( CEREAL_NVP( func_ ) ); // func::FuncOP
	arc( CEREAL_NVP( score_type_ ) ); // const enum core::scoring::ScoreType
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::DOF_Constraint::load( Archive & arc ) {
	arc( const_cast< id::DOF_ID & > ( dof_id_ ) ); // const id::DOF_ID
	arc( func_ ); // func::FuncOP
	arc( const_cast< ScoreType & > ( score_type_ ) ); // const enum core::scoring::ScoreType
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::DOF_Constraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::DOF_Constraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_DOF_Constraint )
#endif // SERIALIZATION
