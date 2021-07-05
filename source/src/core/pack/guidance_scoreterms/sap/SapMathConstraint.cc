// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/guidance_scoreterms/sap/SapMathConstraint.cc
/// @brief A constraint that allows you to subtract and add other SapConstraints
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#include <core/pack/guidance_scoreterms/sap/SapMathConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>

#include <basic/Tracer.hh>


#include <utility/pointer/memory.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR( "core.pack.guidance_scoreterms.sap.SapMathConstraint" );

/// @brief Constructor
///
SapMathConstraint::SapMathConstraint():
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::sap_constraint ),
	upper_bound_( utility::get_undefined_real() ),
	lower_bound_( utility::get_undefined_real() ),
	penalty_per_unit_( 1 )
{}

/// @brief Copy constructor
///
SapMathConstraint::SapMathConstraint( SapMathConstraint const & ot ) :
	core::scoring::aa_composition_energy::SequenceConstraint( core::scoring::sap_constraint )
{
	*this = ot;
}


SapMathConstraint&
SapMathConstraint::operator=( const SapMathConstraint & ot ) {
	helper_weights_ = ot.helper_weights_;
	upper_bound_ = ot.upper_bound_;
	lower_bound_ = ot.lower_bound_;
	penalty_per_unit_ = ot.penalty_per_unit_;
	return *this;
}

/// @brief Destructor
///
SapMathConstraint::~SapMathConstraint() = default;

/// @brief Clone operator
///
core::scoring::constraints::ConstraintOP
SapMathConstraint::clone() const { return utility::pointer::make_shared<SapMathConstraint>( *this ); }

bool SapMathConstraint::operator == ( Constraint const & other ) const
{
	if ( ! other.same_type_as_me( *this ) ) return false;
	if ( !       same_type_as_me( other ) ) return false;

	// Need to actually compare Options, for now just return false
	return false;
}

bool
SapMathConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< SapMathConstraint const * > (&other);
}


void
SapMathConstraint::add_constraint( Real weight, std::string const & name ) {
	helper_weights_.emplace_back( weight, name );
}

void
SapMathConstraint::upper_bound( Real upper ) {
	upper_bound_ = upper;
}

void
SapMathConstraint::lower_bound( Real lower ) {
	lower_bound_ = lower;
}

void
SapMathConstraint::penalty_per_unit( Real penalty ) {
	penalty_per_unit_ = penalty;
}


utility::vector1< std::pair< Real, SapConstraintHelperCOP > >
SapMathConstraint::parse_helpers( utility::vector1< SapConstraintHelperCOP > const & helpers) const {

	utility::vector1< std::pair< Real, SapConstraintHelperCOP > > parsed_helpers;

	for ( std::pair< Real, std::string > weight_name : helper_weights_ ) {

		SapConstraintHelperCOP the_helper = nullptr;
		for ( SapConstraintHelperCOP const & helper : helpers ) {
			if ( helper->options()->name() == weight_name.second ) {
				if ( the_helper ) {
					utility_exit_with_message("SapMathConstraint: Multiple SapConstraints added to pose with name: " + weight_name.second );
				}
				the_helper = helper;
			}
		}
		if ( ! the_helper ) {
			utility_exit_with_message("SapMathConstraint: No SapConstraint added to pose with name: " + weight_name.second );
		}

		parsed_helpers.emplace_back( weight_name.first, the_helper );
	}

	return parsed_helpers;
}


Real
SapMathConstraint::get_math_result( utility::vector1< std::pair< Real, SapConstraintHelperCOP > > const & parsed_helpers ) const {
	Real math_result = 0;

	for ( std::pair< Real, SapConstraintHelperCOP > const & weight_helper : parsed_helpers ) {
		math_result += weight_helper.first * weight_helper.second->current_score();
	}

	return math_result;
}

Real
SapMathConstraint::get_score( utility::vector1< std::pair< Real, SapConstraintHelperCOP > > const & parsed_helpers ) const {

	Real math_result = get_math_result( parsed_helpers );

	if ( ! utility::isnan( lower_bound_ ) && math_result < lower_bound_ ) {
		Real delta = lower_bound_ - math_result;
		return delta * penalty_per_unit_;
	}

	if ( ! utility::isnan( upper_bound_ ) && math_result > upper_bound_ ) {
		Real delta = math_result - upper_bound_;
		return delta * penalty_per_unit_;
	}

	return 0;
}



} //sap
} //guidance_scoreterms
} //pack
} //core

#ifdef    SERIALIZATION
template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapMathConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( CEREAL_NVP( helper_weights_ ) );
	arc( CEREAL_NVP( upper_bound_ ) );
	arc( CEREAL_NVP( lower_bound_ ) );
	arc( CEREAL_NVP( penalty_per_unit_ ) );
}

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapMathConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::aa_composition_energy::SequenceConstraint >( this ) );
	arc( helper_weights_ );
	arc( upper_bound_ );
	arc( lower_bound_ );
	arc( penalty_per_unit_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::guidance_scoreterms::sap::SapMathConstraint );
CEREAL_REGISTER_TYPE( core::pack::guidance_scoreterms::sap::SapMathConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapMathConstraint )
#endif // SERIALIZATION


