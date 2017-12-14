// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/ConstantConstraint.cc

#include <core/types.hh>
#include <core/scoring/constraints/ConstantConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <core/id/AtomID.hh>

#include <utility>
#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @brief Constructor
ConstantConstraint::ConstantConstraint(
	func::FuncOP func_in, // we take ownership of this guy
	ScoreType scotype
):
	Constraint( scotype ),
	func_(std::move( func_in ))
{}

ConstantConstraint::~ConstantConstraint() = default;

ConstraintOP
ConstantConstraint::clone() const {
	// compiler-provided copy constructor performs a shallow copy -- after we
	// use that copy constructor, we need to make a deep copy of the func_
	return ConstantConstraintOP( new ConstantConstraint( *this ));
}

bool ConstantConstraint::operator == ( Constraint const & other ) const {

	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	auto const & other_const( static_cast< ConstantConstraint const & > (other) );

	return func_ == other_const.func_ || ( func_ && other_const.func_ && *func_ == *other_const.func_ );
}

bool ConstantConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< ConstantConstraint const * > (&other);
}

std::string ConstantConstraint::type() const {
	return "Constant";
}

/// @brief compute score
Real
ConstantConstraint::score() const {
	return func_->func(0);
}

/// @brief compute score
void
ConstantConstraint::score( func::XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += score();
}

/// @details In this case, 0 is literally what's being passed to the func
core::Real
ConstantConstraint::dist( core::scoring::func::XYZ_Func const & ) const {
	return 0;
}

/// @brief compute atom deriv
void
ConstantConstraint::fill_f1_f2(
	AtomID const & ,
	func::XYZ_Func const &,
	Vector & /*F1*/,
	Vector & /*F2*/,
	EnergyMap const &
) const
{}

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

ConstantConstraint::ConstantConstraint( ConstantConstraint const & src ) :
	Constraint( src ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::ConstantConstraint::ConstantConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::ConstantConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( func_ ) ); // func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::ConstantConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( func_ ); // func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::ConstantConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::ConstantConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_ConstantConstraint )
#endif // SERIALIZATION
