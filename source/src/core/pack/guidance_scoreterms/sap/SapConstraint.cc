// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/guidance_scoreterms/sap/SapConstraint.cc
/// @brief Apply the sap_score as a sequence constraint to the pose
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#include <core/pack/guidance_scoreterms/sap/SapConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>

#include <basic/Tracer.hh>


#include <utility/pointer/memory.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

#include <utility/exit.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR( "core.pack.guidance_scoreterms.sap.SapConstraint" );

/// @brief Constructor
///
SapConstraint::SapConstraint( SapConstraintOptionsCOP const & options ):
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::sap_constraint ),
	options_( options->clone() )
{}

/// @brief Copy constructor
///
SapConstraint::SapConstraint( SapConstraint const &src ):
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::sap_constraint ),
	options_() //Cloned below.
{
	runtime_assert( src.options_ );
	options_ = src.options_->clone();
}

SapConstraint&
SapConstraint::operator=( const SapConstraint & other ) {
	runtime_assert( other.options_ );
	options_ = other.options_->clone();
	return *this;
}

/// @brief Destructor
///
SapConstraint::~SapConstraint() = default;

/// @brief Clone operator
///
core::scoring::constraints::ConstraintOP
SapConstraint::clone() const { return utility::pointer::make_shared<SapConstraint>( *this ); }

bool SapConstraint::operator == ( Constraint const & other ) const
{
	if ( ! other.same_type_as_me( *this ) ) return false;
	if ( !       same_type_as_me( other ) ) return false;

	// Need to actually compare Options, for now just return false
	return false;
}

bool
SapConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< SapConstraint const * > (&other);
}

SapConstraintOptionsOP
SapConstraint::get_options() {
	return options_;
}

SapConstraintOptionsCOP
SapConstraint::get_const_options() const {
	return options_;
}


} //sap
} //guidance_scoreterms
} //pack
} //core

#ifdef    SERIALIZATION
template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( CEREAL_NVP( options_ ) );
}

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( options_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::guidance_scoreterms::sap::SapConstraint );
CEREAL_REGISTER_TYPE( core::pack::guidance_scoreterms::sap::SapConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapConstraint )
#endif // SERIALIZATION


