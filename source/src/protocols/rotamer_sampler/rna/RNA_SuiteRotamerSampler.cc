// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_SuiteRotamerSampler.cc
/// @brief Generate rotamers for one RNA suite (from residue i to i+1).
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_SuiteRotamerSampler.hh>

// Package headers
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamerSampler.hh>
#include <protocols/rotamer_sampler/RotamerSamplerOneTorsion.hh>
#include <protocols/rotamer_sampler/RotamerSamplerSizedComb.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/rna/util.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_SuiteRotamerSampler" );

namespace protocols {
namespace rotamer_sampler {
namespace rna {

////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_SuiteRotamerSampler::RNA_SuiteRotamerSampler(
	core::Size const rsd_id,
	PuckerState const pucker_state_lower, //ANY_PUCKER, NORTH, SOUTH, NO_PUCKER
	PuckerState const pucker_state_upper,
	ChiState const base_state_lower, //ANY_CHI, ANTI, SYN, NO_CHI
	ChiState const base_state_upper
):
	RotamerSamplerSizedAny(),
	rsd_id_( rsd_id ),
	base_state_lower_( base_state_lower ),
	base_state_upper_( base_state_upper ),
	sample_nucleoside_lower_( true ),
	sample_nucleoside_upper_( true ),
	extra_epsilon_( false ),
	extra_beta_( false ),
	extra_chi_( false ),
	skip_same_pucker_( true ),
	idealize_coord_( true ),
	fast_( false ),
	bin_size_( 20 )
{
	runtime_assert( pucker_state_lower <= 2 );
	runtime_assert( pucker_state_upper <= 2 );
	runtime_assert( base_state_lower <= 3 );
	runtime_assert( base_state_upper <= 3 );

	if ( pucker_state_lower == ANY_PUCKER ) {
		pucker_states_lower_.push_back( NORTH );
		pucker_states_lower_.push_back( SOUTH );
	} else {
		pucker_states_lower_.push_back( pucker_state_lower );
	}

	if ( pucker_state_upper == ANY_PUCKER ) {
		pucker_states_upper_.push_back( NORTH );
		pucker_states_upper_.push_back( SOUTH );
	} else {
		pucker_states_upper_.push_back( pucker_state_upper );
	}
}
////////////////////////////////////////////////////////////////////////////
void RNA_SuiteRotamerSampler::init() {
	runtime_assert( pucker_states_lower_.size() == 1 ||
			sample_nucleoside_lower_ );
	runtime_assert( pucker_states_upper_.size() == 1 ||
			sample_nucleoside_upper_ );
	clear_rotamer();

	if ( fast_ ) {
		init_fast();
	} else {
		init_standard();
	}

	RotamerSamplerSizedAny::init();
}
//////////////////////////////////////////////////////////////////////////
void RNA_SuiteRotamerSampler::init_standard() {
	using namespace core::id;
	TorsionList full_torsions = get_full_torsions( bin_size_ );

	RNA_FittedTorsionInfo const torsion_info;
	//Setup the rotamer samplers
	for ( Size i = 1; i <= pucker_states_lower_.size(); ++i ) {
		for ( Size j = 1; j <= pucker_states_upper_.size(); ++j ) {
			//Join individual torsions combinatorially.
			RotamerSamplerSizedCombOP new_rotamer_agg = new RotamerSamplerSizedComb;

			/////Gamma rotamers///// (full torsion)
			RotamerSamplerOneTorsionOP gamma_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, GAMMA ), full_torsions );
			new_rotamer_agg->add_external_loop_rotamer( gamma_rotamer );

			/////Beta rotamers/////
			//default: 180 +- 100 deg
			//extra_beta: sample full torsion
			TorsionList beta_torsions;
			if ( extra_beta_ ) {
				beta_torsions = full_torsions;
			} else {
				add_values_from_center( beta_torsions, 180, 100, bin_size_ );
			}
			RotamerSamplerOneTorsionOP beta_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, BETA ), beta_torsions );
			new_rotamer_agg->add_external_loop_rotamer( beta_rotamer );

			/////Alpha rotamers///// (full torsion)
			RotamerSamplerOneTorsionOP alpha_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, ALPHA ), full_torsions );
			new_rotamer_agg->add_external_loop_rotamer( alpha_rotamer );

			/////Zeta rotamers///// (full torsion)
			RotamerSamplerOneTorsionOP zeta_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_, BB, ZETA ), full_torsions );
			new_rotamer_agg->add_external_loop_rotamer( zeta_rotamer );

			/////Epsilon rotamers/////
			//default: center +- 20 deg
			//extra_epsilon: center +- 60 deg
			TorsionList epsilon_torsions = get_epsilon_torsions( (pucker_states_lower_[i] == NORTH), extra_epsilon_, bin_size_ );
			RotamerSamplerOneTorsionOP epsilon_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_, BB, EPSILON ), epsilon_torsions );
			new_rotamer_agg->add_external_loop_rotamer( epsilon_rotamer );

			/////Nucleoside rotamers/////
			if ( sample_nucleoside_lower_ ) {
				RNA_NucleosideRotamerSamplerOP nucleoside_rotamer1 =
						new RNA_NucleosideRotamerSampler( rsd_id_,
						pucker_states_lower_[i], base_state_lower_ );
				nucleoside_rotamer1->set_bin_size( bin_size_ );
				nucleoside_rotamer1->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer1->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer1->set_extra_chi( extra_chi_ );
				new_rotamer_agg->add_external_loop_rotamer( nucleoside_rotamer1 );
			}

			if ( sample_nucleoside_upper_ ) {
				RNA_NucleosideRotamerSamplerOP nucleoside_rotamer2 =
						new RNA_NucleosideRotamerSampler( rsd_id_ + 1,
						pucker_states_upper_[j], base_state_upper_ );
				nucleoside_rotamer2->set_bin_size( bin_size_ );
				nucleoside_rotamer2->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer2->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer2->set_extra_chi( extra_chi_ );
				new_rotamer_agg->add_external_loop_rotamer( nucleoside_rotamer2 );
			}

			/////Add to this sampler/////
			add_external_loop_rotamer( new_rotamer_agg );
		}
	}
}
//////////////////////////////////////////////////////////////////////////
//Fast sampling for integration test, only sample center torsions
void RNA_SuiteRotamerSampler::init_fast() {
	using namespace core::id;
	extra_epsilon_ = false;
	extra_beta_ = false;
	extra_chi_ = false;

	RNA_FittedTorsionInfo const torsion_info;
	TorsionList torsions;

	//Setup the rotamer samplers
	for ( Size i = 1; i <= pucker_states_lower_.size(); ++i ) {
		for ( Size j = 1; j <= pucker_states_upper_.size(); ++j ) {
			//Join individual torsions combinatorially.
			RotamerSamplerSizedCombOP new_rotamer_agg = new RotamerSamplerSizedComb;

			/////Gamma rotamers/////
			torsions = fast_sample_torsions_from_info(
					torsion_info.gaussian_parameter_set_gamma() );
			RotamerSamplerOneTorsionOP gamma_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, GAMMA ), torsions );
			new_rotamer_agg->add_external_loop_rotamer( gamma_rotamer );

			/////Beta rotamers/////
			torsions = fast_sample_torsions_from_info(
					torsion_info.gaussian_parameter_set_beta() );
			RotamerSamplerOneTorsionOP beta_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, BETA ), torsions );
			new_rotamer_agg->add_external_loop_rotamer( beta_rotamer );

			/////Alpha rotamers/////
			torsions = fast_sample_torsions_from_info(
					torsion_info.gaussian_parameter_set_alpha() );
			RotamerSamplerOneTorsionOP alpha_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, ALPHA ), torsions );
			new_rotamer_agg->add_external_loop_rotamer( alpha_rotamer );

			/////Zeta rotamers/////
			torsions = fast_sample_torsions_from_info(
					torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus() );
			RotamerSamplerOneTorsionOP zeta_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_, BB, ZETA ), torsions );
			new_rotamer_agg->add_external_loop_rotamer( zeta_rotamer );

			/////Epsilon rotamers/////
			if ( pucker_states_lower_[i] == NORTH ) {
				torsions = fast_sample_torsions_from_info(
						torsion_info.gaussian_parameter_set_epsilon_north() );
			} else {
				torsions = fast_sample_torsions_from_info(
						torsion_info.gaussian_parameter_set_epsilon_south() );
			}
			RotamerSamplerOneTorsionOP epsilon_rotamer = new RotamerSamplerOneTorsion(
					TorsionID( rsd_id_, BB, EPSILON ), torsions );
			new_rotamer_agg->add_external_loop_rotamer( epsilon_rotamer );

			/////Nucleoside rotamers/////
			if ( sample_nucleoside_lower_ ) {
				RNA_NucleosideRotamerSamplerOP nucleoside_rotamer1 =
						new RNA_NucleosideRotamerSampler( rsd_id_,
						pucker_states_lower_[i], base_state_lower_ );
				nucleoside_rotamer1->set_fast( fast_ );
				nucleoside_rotamer1->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer1->set_idealize_coord( idealize_coord_ );
				new_rotamer_agg->add_external_loop_rotamer( nucleoside_rotamer1 );
			}

			if ( sample_nucleoside_upper_ ) {
				RNA_NucleosideRotamerSamplerOP nucleoside_rotamer2 =
						new RNA_NucleosideRotamerSampler( rsd_id_ + 1,
						pucker_states_upper_[j], base_state_upper_ );
				nucleoside_rotamer2->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer2->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer2->set_fast( fast_ );
				new_rotamer_agg->add_external_loop_rotamer( nucleoside_rotamer2 );
			}

			/////Add to this sampler/////
			add_external_loop_rotamer( new_rotamer_agg );
		}
	}
}
//////////////////////////////////////////////////////////////////////////
RNA_SuiteRotamerSampler::TorsionList
RNA_SuiteRotamerSampler::fast_sample_torsions_from_info(
	utility::vector1<GaussianParameter> const & params
) {
	TorsionList torsions;
	for ( Size i = 1; i <= params.size(); ++i ) {
		for ( int k = -1; k <= 1; k++ ){
			torsions.push_back( params[i].center + k * params[i].width );
		}
	}
	return torsions;
}
//////////////////////////////////////////////////////////////////////////

} //rna
} //rotamer_sampler
} //protocols
