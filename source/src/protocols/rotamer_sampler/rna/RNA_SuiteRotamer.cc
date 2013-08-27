// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_SuiteRotamer.cc
/// @brief Generate rotamers for one RNA suite (from residue i to i+1).
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_SuiteRotamer.hh>

// Package headers
#include <protocols/rotamer_sampler/rna/RNA_NucleosideRotamer.hh>
#include <protocols/rotamer_sampler/RotamerOneTorsion.hh>
#include <protocols/rotamer_sampler/RotamerSizedComb.hh>
#include <protocols/rotamer_sampler/util.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_SuiteRotamer" );

namespace protocols {
namespace rotamer_sampler {
namespace rna {

////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_SuiteRotamer::RNA_SuiteRotamer(
	core::Size const rsd_id,
	core::Size const pucker_state_lower, //WHATEVER, NORTH, SOUTH
	core::Size const pucker_state_upper,
	core::Size const base_state_lower, //WHATEVER, ANTI, SYN, NONE
	core::Size const base_state_upper
):
	RotamerSizedAny(),
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

	if ( pucker_state_lower == WHATEVER ) {
		pucker_states_lower_.push_back( NORTH );
		pucker_states_lower_.push_back( SOUTH );
	} else {
		pucker_states_lower_.push_back( pucker_state_lower );
	}

	if ( pucker_state_upper == WHATEVER ) {
		pucker_states_upper_.push_back( NORTH );
		pucker_states_upper_.push_back( SOUTH );
	} else {
		pucker_states_upper_.push_back( pucker_state_upper );
	}
}
////////////////////////////////////////////////////////////////////////////
void RNA_SuiteRotamer::init() {
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

	RotamerSizedAny::init();
}
//////////////////////////////////////////////////////////////////////////
void RNA_SuiteRotamer::init_standard() {
	using namespace core::id;
	TorsionList full_torsions;
	add_values_from_center( full_torsions, 0, 180, bin_size_ );
	//Avoid sampling both -180 and 180 deg
	Real const epsil =  0.001; //Arbitary small number
	if ( 360 - ( full_torsions.back() - full_torsions.front() ) < epsil )
			full_torsions.pop_back();

	RNA_FittedTorsionInfo const torsion_info;
	//Setup the rotamer samplers
	for ( Size i = 1; i <= pucker_states_lower_.size(); ++i ) {
		for ( Size j = 1; j <= pucker_states_upper_.size(); ++j ) {
			//Join individual torsions combinatorially.
			RotamerSizedCombOP new_rotamer_agg = new RotamerSizedComb;

			/////Gamma rotamers///// (full torsion)
			RotamerOneTorsionOP gamma_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, GAMMA ), full_torsions );
			new_rotamer_agg->add_rotamer( gamma_rotamer );

			/////Beta rotamers/////
			//default: 180 +- 100 deg
			//extra_beta: sample full torsion
			TorsionList beta_torsions;
			if ( extra_beta_ ) {
				beta_torsions = full_torsions;
			} else {
				add_values_from_center( beta_torsions, 180, 100, bin_size_ );
			}
			RotamerOneTorsionOP beta_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, BETA ), beta_torsions );
			new_rotamer_agg->add_rotamer( beta_rotamer );

			/////Alpha rotamers///// (full torsion)
			RotamerOneTorsionOP alpha_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, ALPHA ), full_torsions );
			new_rotamer_agg->add_rotamer( alpha_rotamer );

			/////Zeta rotamers///// (full torsion)
			RotamerOneTorsionOP zeta_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_, BB, ZETA ), full_torsions );
			new_rotamer_agg->add_rotamer( zeta_rotamer );

			/////Epsilon rotamers/////
			//default: center +- 20 deg
			//extra_epsilon: center +- 60 deg
			TorsionList epsilon_torsions;
			Real center = ( pucker_states_lower_[i] == NORTH ) ?
					torsion_info.epsilon_north() : torsion_info.epsilon_south();
			Real max_range = 20;
			if ( extra_epsilon_ ) {
				max_range = 60;
				//Choice made by Parin, to cover the uneven tails of
				//epsilons distributions by extra_epsilon mode.
				if ( pucker_states_lower_[i] == NORTH ) {
					center -= 20;
				} else {
					center += 20;
				}
			}
			add_values_from_center( epsilon_torsions, center, max_range, bin_size_ );
			RotamerOneTorsionOP epsilon_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_, BB, EPSILON ), epsilon_torsions );
			new_rotamer_agg->add_rotamer( epsilon_rotamer );

			/////Nucleoside rotamers/////
			if ( sample_nucleoside_lower_ ) {
				RNA_NucleosideRotamerOP nucleoside_rotamer1 =
						new RNA_NucleosideRotamer( rsd_id_,
						pucker_states_lower_[i], base_state_lower_ );
				nucleoside_rotamer1->set_bin_size( bin_size_ );
				nucleoside_rotamer1->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer1->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer1->set_extra_chi( extra_chi_ );
				new_rotamer_agg->add_rotamer( nucleoside_rotamer1 );
			}

			if ( sample_nucleoside_upper_ ) {
				RNA_NucleosideRotamerOP nucleoside_rotamer2 =
						new RNA_NucleosideRotamer( rsd_id_ + 1,
						pucker_states_upper_[j], base_state_upper_ );
				nucleoside_rotamer2->set_bin_size( bin_size_ );
				nucleoside_rotamer2->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer2->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer2->set_extra_chi( extra_chi_ );
				new_rotamer_agg->add_rotamer( nucleoside_rotamer2 );
			}

			/////Add to this sampler/////
			add_rotamer( new_rotamer_agg );
		}
	}
}
//////////////////////////////////////////////////////////////////////////
//Fast sampling for integration test, only sample center torsions
void RNA_SuiteRotamer::init_fast() {
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
			RotamerSizedCombOP new_rotamer_agg = new RotamerSizedComb;

			/////Gamma rotamers/////
			torsions = torsions_from_info(
					torsion_info.gaussian_parameter_set_gamma() );
			RotamerOneTorsionOP gamma_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, GAMMA ), torsions );
			new_rotamer_agg->add_rotamer( gamma_rotamer );

			/////Beta rotamers/////
			torsions = torsions_from_info(
					torsion_info.gaussian_parameter_set_beta() );
			RotamerOneTorsionOP beta_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, BETA ), torsions );
			new_rotamer_agg->add_rotamer( beta_rotamer );

			/////Alpha rotamers/////
			torsions = torsions_from_info(
					torsion_info.gaussian_parameter_set_alpha() );
			RotamerOneTorsionOP alpha_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_ + 1, BB, ALPHA ), torsions );
			new_rotamer_agg->add_rotamer( alpha_rotamer );

			/////Zeta rotamers/////
			torsions = torsions_from_info(
					torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus() );
			RotamerOneTorsionOP zeta_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_, BB, ZETA ), torsions );
			new_rotamer_agg->add_rotamer( zeta_rotamer );

			/////Epsilon rotamers/////
			if ( pucker_states_lower_[i] == NORTH ) {
				torsions = torsions_from_info(
						torsion_info.gaussian_parameter_set_epsilon_north() );
			} else {
				torsions = torsions_from_info(
						torsion_info.gaussian_parameter_set_epsilon_south() );
			}
			RotamerOneTorsionOP epsilon_rotamer = new RotamerOneTorsion(
					TorsionID( rsd_id_, BB, EPSILON ), torsions );
			new_rotamer_agg->add_rotamer( epsilon_rotamer );

			/////Nucleoside rotamers/////
			if ( sample_nucleoside_lower_ ) {
				RNA_NucleosideRotamerOP nucleoside_rotamer1 =
						new RNA_NucleosideRotamer( rsd_id_,
						pucker_states_lower_[i], base_state_lower_ );
				nucleoside_rotamer1->set_fast( fast_ );
				nucleoside_rotamer1->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer1->set_idealize_coord( idealize_coord_ );
				new_rotamer_agg->add_rotamer( nucleoside_rotamer1 );
			}

			if ( sample_nucleoside_upper_ ) {
				RNA_NucleosideRotamerOP nucleoside_rotamer2 =
						new RNA_NucleosideRotamer( rsd_id_ + 1,
						pucker_states_upper_[j], base_state_upper_ );
				nucleoside_rotamer2->set_skip_same_pucker( skip_same_pucker_ );
				nucleoside_rotamer2->set_idealize_coord( idealize_coord_ );
				nucleoside_rotamer2->set_fast( fast_ );
				new_rotamer_agg->add_rotamer( nucleoside_rotamer2 );
			}

			/////Add to this sampler/////
			add_rotamer( new_rotamer_agg );
		}
	}
}
//////////////////////////////////////////////////////////////////////////
RNA_SuiteRotamer::TorsionList
RNA_SuiteRotamer::torsions_from_info(
	utility::vector1<GaussianParameter> const & params
) {
	TorsionList torsions;
	for ( Size i = 1; i <= params.size(); ++i ) {
		torsions.push_back( params[i].center );
	}
	return torsions;
}
//////////////////////////////////////////////////////////////////////////

}
}
}
