// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.cc
/// @brief Generate glycosidic chi rotamers for RNA.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>

// Package headers

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <basic/Tracer.hh>


using namespace core;
using namespace core::chemical::rna;
static basic::Tracer TR( "protocols.sampler.rna.RNA_ChiStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

//////////////////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_ChiStepWiseSampler::RNA_ChiStepWiseSampler(
	Size const rsd_id,
	PuckerState const pucker_state,
	ChiState const base_state
):
	StepWiseSamplerOneTorsion(),
	rsd_id_( rsd_id ),
	base_state_( base_state ),
	pucker_state_( pucker_state ),
	bin_size_( 20 ),
	max_range_( 20 ) //+- 20 degrees
{
	runtime_assert( base_state <= 2 );
	runtime_assert( pucker_state == NORTH || pucker_state_ == SOUTH );
}
//////////////////////////////////////////////////////////////////////////
void RNA_ChiStepWiseSampler::init() {
	Real chi_center;
	TorsionList allowed_torions;
	if ( pucker_state_ == NORTH ) {
		if ( base_state_ == ANTI || base_state_ == ANY_CHI ) {
			chi_center = torsion_info_.chi_north_anti();
			add_values_from_center( allowed_torions, chi_center, max_range_, bin_size_ );
		}
		if ( base_state_ == SYN || base_state_ == ANY_CHI ) {
			chi_center = torsion_info_.chi_north_syn();
			add_values_from_center( allowed_torions, chi_center, max_range_, bin_size_ );
		}
	} else if ( pucker_state_ == SOUTH ) {
		if ( base_state_ == ANTI || base_state_ == ANY_CHI ) {
			chi_center = torsion_info_.chi_south_anti();
			add_values_from_center( allowed_torions, chi_center, max_range_, bin_size_ );
		}
		if ( base_state_ == SYN || base_state_ == ANY_CHI ) {
			chi_center = torsion_info_.chi_south_syn();
			add_values_from_center( allowed_torions, chi_center, max_range_, bin_size_ );
		}
	}

	set_torsion_id( id::TorsionID( rsd_id_, id::CHI, 1 ) );
	set_torsions( allowed_torions );

	StepWiseSamplerOneTorsion::init();
}
} //rna
} //sampler
} //stepwise
} //protocols
