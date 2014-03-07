// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_McSugar.cc
/// @brief Markov chain sampler for sugar pucker.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_McSugar.hh>

// Package headers
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
using namespace core::chemical::rna;

static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_McSugar" );
static numeric::random::RandomGenerator RG( 4864698 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
namespace rna {

///////////////////////////////////////////////////////////////////////////
RNA_McSugar::RNA_McSugar(
	Size const rsd_id,
	Real const flip_rate,
	PuckerState const init_pucker
):
	McRotamer(),
	stored_pucker_state_( init_pucker ),
	active_pucker_state_( init_pucker ),
	flip_rate_( flip_rate ),
	sugar_rotamer_( new RNA_SugarRotamer( rsd_id, ANY_PUCKER ) )
{}
///////////////////////////////////////////////////////////////////////////
void RNA_McSugar::init() {
	set_init( true );
	sugar_rotamer_->set_random( true );
	sugar_rotamer_->init();
	if ( uniform_sampling() ) flip_rate_ = 1;
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RNA_McSugar::operator++() {
	runtime_assert( is_init() );
	if ( flip_rate_ != 0 && RG.uniform() < flip_rate_ ) {
		active_pucker_state_ = ( stored_pucker_state_ == NORTH ) ? SOUTH : NORTH;
	} else {
		active_pucker_state_ = stored_pucker_state_;
	}
}
///////////////////////////////////////////////////////////////////////////
void RNA_McSugar::apply( pose::Pose & pose ) {
	runtime_assert( is_init() );
	sugar_rotamer_->apply( pose, active_pucker_state_ );
}
///////////////////////////////////////////////////////////////////////////
void RNA_McSugar::set_skip_same_pucker( bool const setting ) {
	set_init( false );
	sugar_rotamer_->set_skip_same_pucker( setting );
}
///////////////////////////////////////////////////////////////////////////
void RNA_McSugar::set_idealize_coord( bool const setting ) {
	set_init( false );
	sugar_rotamer_->set_idealize_coord( setting );
}
///////////////////////////////////////////////////////////////////////////
}
}
}
