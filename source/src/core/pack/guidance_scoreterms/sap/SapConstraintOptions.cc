// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh
/// @brief Options for SapConstraint
/// @details Contains all the options and a method to transform sap_score into sap_constraint
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

// Options system

// File I/O

// Other Headers
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>

// C++ headers

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR("core.pack.guidance_scoreterms.sap.SapConstraintOptions");


SapConstraintOptions::SapConstraintOptions() :
	utility::VirtualBase(),
	sap_goal_( 0 ),
	sap_lb_goal_( 0 ),
	packing_correction_( 0 ),
	penalty_per_sap_( 1 ),
	score_selector_( utility::pointer::make_shared< select::residue_selector::TrueResidueSelector>() ),
	sap_calculate_selector_( nullptr ),
	sasa_selector_( nullptr ),
	fast_( false ),
	lightning_( false ),
	full_accuracy_when_scoring_( true ),
	name_( "" )
{}

SapConstraintOptionsOP SapConstraintOptions::clone() const {
	return utility::pointer::make_shared<SapConstraintOptions>( *this );
}


void
SapConstraintOptions::sap_goal( Real goal ) {
	sap_goal_ = goal;
}

void
SapConstraintOptions::sap_lb_goal( Real lb_goal ) {
	sap_lb_goal_ = lb_goal;
}

void
SapConstraintOptions::packing_correction( Real adjust ) {
	packing_correction_ = adjust;
}

void
SapConstraintOptions::penalty_per_sap( Real penalty ) {
	penalty_per_sap_ = penalty;
}

void
SapConstraintOptions::score_selector( select::residue_selector::ResidueSelectorCOP const & sel ) {
	if ( ! sel ) {
		utility_exit_with_message("SapConstraintOptions: You can't pass a nullptr to score_selector(). It defaults to TrueSelector anyways.");
	}
	score_selector_ = sel->clone();
}

void
SapConstraintOptions::sap_calculate_selector( select::residue_selector::ResidueSelectorCOP const & sel ) {
	if ( ! sel ) {
		sap_calculate_selector_ = nullptr;
	}
	sap_calculate_selector_ = sel->clone();
}

void
SapConstraintOptions::sasa_selector( select::residue_selector::ResidueSelectorCOP const & sel ) {
	if ( ! sel ) {
		sasa_selector_ = nullptr;
	}
	sasa_selector_ = sel->clone();
}

void
SapConstraintOptions::fast( bool fast ) {
	fast_ = fast;
}

void
SapConstraintOptions::lightning( bool lightning ) {
	lightning_ = lightning;
}

void
SapConstraintOptions::name( std::string const & name ) {
	name_ = name;
}


void
SapConstraintOptions::sanity_check( ) const {
	if ( sap_goal() < sap_lb_goal() ) {
		utility_exit_with_message("sap_goal must be greater than or equal to sap_lb_goal!");
	}
}

void
SapConstraintOptions::full_accuracy_when_scoring( bool full_accuracy ) {
	full_accuracy_when_scoring_ = full_accuracy;
}

bool
SapConstraintOptions::full_accuracy_when_scoring() const {
	return full_accuracy_when_scoring_;
}

Real
SapConstraintOptions::sap_goal() const {
	return sap_goal_;
}

Real
SapConstraintOptions::sap_lb_goal() const {
	return sap_lb_goal_;
}

Real
SapConstraintOptions::packing_correction() const {
	return packing_correction_;
}

Real
SapConstraintOptions::penalty_per_sap() const {
	return penalty_per_sap_;
}

std::string
SapConstraintOptions::name() const {
	return name_;
}

select::residue_selector::ResidueSelectorCOP
SapConstraintOptions::score_selector() const {
	return score_selector_;
}

select::residue_selector::ResidueSelectorCOP
SapConstraintOptions::sap_calculate_selector() const {
	if ( sap_calculate_selector_ ) return sap_calculate_selector_;
	return score_selector_;
}

select::residue_selector::ResidueSelectorCOP
SapConstraintOptions::sasa_selector() const {
	if ( sasa_selector_ ) return sasa_selector_;
	if ( sap_calculate_selector_ ) return sap_calculate_selector_;
	return score_selector_;
}

bool
SapConstraintOptions::fast() const {
	return fast_ || lightning_;
}

bool
SapConstraintOptions::lightning() const {
	return lightning_;
}

Real
SapConstraintOptions::transform_sap_to_score( Real score, bool /**/ ) const {
	score += packing_correction();
	if ( score >= sap_goal() ) {
		return ( score - sap_goal() ) * penalty_per_sap();
	}
	if ( score <= sap_lb_goal() ) {
		return ( sap_lb_goal() - score ) * penalty_per_sap();
	}
	return 0;

}



} //sap
} //guidance_scoreterms
} //pack
} //core


#ifdef SERIALIZATION

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapConstraintOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( sap_goal_ ) );
	arc( CEREAL_NVP( sap_lb_goal_ ) );
	arc( CEREAL_NVP( packing_correction_ ) );
	arc( CEREAL_NVP( penalty_per_sap_ ) );
	arc( CEREAL_NVP( score_selector_ ) );
	arc( CEREAL_NVP( sap_calculate_selector_ ) );
	arc( CEREAL_NVP( sasa_selector_ ) );
	arc( CEREAL_NVP( fast_ ) );
	arc( CEREAL_NVP( lightning_ ) );
	arc( CEREAL_NVP( full_accuracy_when_scoring_ ) );
	arc( CEREAL_NVP( name_ ) );
}

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapConstraintOptions::load( Archive & arc ) {
	arc( sap_goal_ );
	arc( sap_lb_goal_ );
	arc( packing_correction_ );
	arc( penalty_per_sap_ );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector );
	score_selector_ = local_selector;

	arc( local_selector );
	sap_calculate_selector_ = local_selector;

	arc( local_selector );
	sasa_selector_ = local_selector;

	arc( fast_ );
	arc( lightning_ );
	arc( full_accuracy_when_scoring_ );
	arc( name_ );
}


SAVE_AND_LOAD_SERIALIZABLE( core::pack::guidance_scoreterms::sap::SapConstraintOptions );
CEREAL_REGISTER_TYPE( core::pack::guidance_scoreterms::sap::SapConstraintOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapConstraintOptions )
#endif // SERIALIZATION
