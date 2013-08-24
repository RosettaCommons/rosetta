// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_SugarRotamer.cc
/// @brief Generate sugar pucker rotamers for RNA.
/// @detailed
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.hh>

// Package headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_SugarRotamer" );
static numeric::random::RandomGenerator RG( 256199 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_SugarRotamer::RNA_SugarRotamer(
	Size const rsd_id,
	Size const pucker_state
):
	RotamerSized(),
	id_( 0 ),
	rsd_id_( rsd_id ),
	pucker_state_( pucker_state ),
	skip_same_pucker_( true ),
	idealize_coord_( true )
{
	runtime_assert( pucker_state_ == WHATEVER || pucker_state_ == NORTH || pucker_state_ == SOUTH );
}
//////////////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::init() {
	set_init( true );
	if ( pucker_state_ == WHATEVER ) {
		pucker_states_.push_back( NORTH );
		pucker_states_.push_back( SOUTH );
	} else {
		pucker_states_.push_back( pucker_state_ );
	}
	reset();
}
//////////////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::reset() {
	runtime_assert( is_init() );
	if ( is_random() ) {
		++( *this );
	} else {
		id_ = 1;
	}
}
//////////////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::operator++() {
	runtime_assert( not_end() );
	if ( is_random() ) {
		id_ = RG.random_range( 1, size() );
	} else {
		++id_;
	}
}
///////////////////////////////////////////////////////////////////////////
bool RNA_SugarRotamer::not_end() const {
	runtime_assert( is_init() );
	if ( is_random() ) return true;
	return ( id_ <= size() );
}
//////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::apply( pose::Pose & pose, core::Size const i ) {
	runtime_assert( is_init() );
	apply_pucker(	pose, rsd_id_, pucker_states_[i],
		skip_same_pucker_, idealize_coord_ );
}
//////////////////////////////////////////////////////////////////
}
}
}
