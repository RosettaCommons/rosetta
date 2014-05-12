// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.cc
/// @brief Generate rotamers for one RNA nucleoside (pucker + glycosidic chi).
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.hh>

// Package headers
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.hh>
#include <protocols/rotamer_sampler/rna/RNA_ChiRotamer.hh>
#include <protocols/rotamer_sampler/RotamerSizedComb.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/util.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_NucleosideRotamer" );

namespace protocols {
namespace rotamer_sampler {
namespace rna {

////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_NucleosideRotamer::RNA_NucleosideRotamer(
	core::Size const rsd_id,
	PuckerState const pucker_state, //ANY_PUCKER, NORTH, SOUTH, NO_PUCKER
	ChiState const base_state //ANY_CHI, ANTI, SYN, NO_CHI
):
	RotamerSizedAny(),
	rsd_id_( rsd_id ),
	base_state_( base_state ),
	extra_chi_( false ),
	skip_same_pucker_( true ),
	idealize_coord_( true ),
	fast_( false ),
	bin_size_( 20 )
{
	runtime_assert( pucker_state <= 2 );
	runtime_assert( base_state <= 3 );

	if ( pucker_state == ANY_PUCKER ) {
		pucker_states_.push_back( NORTH );
		pucker_states_.push_back( SOUTH );
	} else {
		pucker_states_.push_back( pucker_state );
	}
}
////////////////////////////////////////////////////////////////////////////
void RNA_NucleosideRotamer::init() {
	clear_rotamer();

	//Setup the rotamer samplers
	for ( Size i = 1; i <= pucker_states_.size(); ++i ) {
		RotamerSizedCombOP new_rotamer_agg = new RotamerSizedComb;

		/////Chi rotamers/////
		if ( base_state_ != NO_CHI ) {
			RNA_ChiRotamerOP chi_rotamer = new RNA_ChiRotamer(
					rsd_id_, pucker_states_[i], base_state_ );
			chi_rotamer->set_bin_size( bin_size_ );
			chi_rotamer->set_extra_chi( extra_chi_);
			if ( fast_ ) chi_rotamer->set_max_range( 1 );
			new_rotamer_agg->add_rotamer( chi_rotamer );
		}

		/////Pucker rotamers/////
		RNA_SugarRotamerOP pucker_rotamer = new RNA_SugarRotamer(
				rsd_id_, pucker_states_[i] );
		pucker_rotamer->set_skip_same_pucker( skip_same_pucker_ );
		pucker_rotamer->set_idealize_coord( idealize_coord_ );
		new_rotamer_agg->add_rotamer( pucker_rotamer );

		/////Add to the this sampler/////
		add_rotamer( new_rotamer_agg );
	}

	RotamerSizedAny::init();
}
//////////////////////////////////////////////////////////////////////////
}
}
}
