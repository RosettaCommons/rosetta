// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_MC_Suite.cc
/// @brief Markov chain sampler for RNA suite.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>

// Package headers
#include <protocols/stepwise/sampler/rna/RNA_MC_Sugar.hh>
#include <protocols/stepwise/sampler/MC_OneTorsion.hh>

// Project headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;

static thread_local basic::Tracer TR( "protocols.sampler.rna.RNA_MC_Suite" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

///////////////////////////////////////////////////////////////////////////
RNA_MC_Suite::RNA_MC_Suite( Size const rsd_id ):
	MC_Comb(),
	rsd_id_( rsd_id ),
	skip_same_pucker_( true ),
	idealize_coord_( true ),
	sample_near_a_form_( false ),
	sample_bb_( true ),
	sample_lower_nucleoside_( false ),
	sample_upper_nucleoside_( true ),
	pucker_flip_rate_( 0.1 ),
	gaussian_stdev_( 20 ),
	a_form_range_( 60 ),
	init_pucker_( NORTH )
{
	using namespace core::id;

	RNA_FittedTorsionInfo const torsion_info;
	a_form_torsions_.clear();
	a_form_torsions_.push_back( torsion_info.gamma_aform() );
	a_form_torsions_.push_back( torsion_info.beta_aform() );
	a_form_torsions_.push_back( torsion_info.alpha_aform() );
	a_form_torsions_.push_back( torsion_info.epsilon_aform() );
	a_form_torsions_.push_back( torsion_info.zeta_aform() );
	a_form_torsions_.push_back( torsion_info.chi_north_anti() );
	a_form_torsions_.push_back( torsion_info.chi_north_anti() );
	init_torsions_ = a_form_torsions_;

	torsion_ids_.clear();
	torsion_ids_.push_back( TorsionID( rsd_id + 1, BB, GAMMA ) );
	torsion_ids_.push_back( TorsionID( rsd_id + 1, BB, BETA ) );
	torsion_ids_.push_back( TorsionID( rsd_id + 1, BB, ALPHA ) );
	torsion_ids_.push_back( TorsionID( rsd_id, BB, EPSILON ) );
	torsion_ids_.push_back( TorsionID( rsd_id, BB, ZETA ) );
	torsion_ids_.push_back( TorsionID( rsd_id, id::CHI, 1 ) );
	torsion_ids_.push_back( TorsionID( rsd_id + 1, id::CHI, 1 ) );
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_Suite::init() {
	set_init( true );
	clear_rotamer();
	if ( sample_near_a_form_ ) {
		pucker_flip_rate_ = 0;
		init_pucker_ = NORTH;
	}

	for ( Size i = 1; i <= 5; ++i ) {
		MC_OneTorsionOP sampler( new MC_OneTorsion( torsion_ids_[i], init_torsions_[i] ) );
		sampler->set_gaussian_stdev( gaussian_stdev_ );
		if ( sample_near_a_form_ ) {
			Real const min_angle( a_form_torsions_[i] - a_form_range_ );
			Real const max_angle( a_form_torsions_[i] + a_form_range_ );
			sampler->set_angle_range( min_angle, max_angle );
		}
		bb_samplers_.push_back( sampler );
	}

	for ( Size i = 1; i <= 2; ++i ) {
		MC_OneTorsionOP chi_sampler( new MC_OneTorsion(
			torsion_ids_[5 + i], init_torsions_[5 + i] ) );
		chi_sampler->set_gaussian_stdev( gaussian_stdev_ );
		if ( sample_near_a_form_ ) {
			Real const min_angle( a_form_torsions_[5 + i] - a_form_range_ );
			Real const max_angle( a_form_torsions_[5 + i] + a_form_range_ );
			chi_sampler->set_angle_range( min_angle, max_angle );
			chi_sampler->set_angle_range( min_angle, max_angle );
		}
		RNA_MC_SugarOP sugar_sampler( new RNA_MC_Sugar(
			rsd_id_ - 1 + i, pucker_flip_rate_, init_pucker_ ) );
		sugar_sampler->set_skip_same_pucker( skip_same_pucker_ );
		sugar_sampler->set_idealize_coord( idealize_coord_ );
		chi_samplers_.push_back( chi_sampler );
		sugar_samplers_.push_back( sugar_sampler );
	}

	if ( sample_bb_ ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++ i ) {
			add_external_loop_rotamer( bb_samplers_[i] );
		}
	}

	if ( sample_lower_nucleoside_ ) {
		add_external_loop_rotamer( chi_samplers_[1] );
		add_external_loop_rotamer( sugar_samplers_[1] );
	}

	if ( sample_upper_nucleoside_ ) {
		add_external_loop_rotamer( chi_samplers_[2] );
		add_external_loop_rotamer( sugar_samplers_[2] );
	}

	MC_Comb::init();
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_Suite::set_init_from_pose( pose::Pose const & pose ) {
	for ( Size i = 1; i <= init_torsions_.size(); ++i ) {
		init_torsions_[i] = pose.torsion( torsion_ids_[i] );
	}
	init_pucker_ = assign_pucker( pose, rsd_id_ );
	set_init( false );
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_Suite::set_pucker_flip_rate( Real const setting ) {
	pucker_flip_rate_ = setting;
	if ( is_init() ) {
		sugar_samplers_[1]->set_flip_rate( pucker_flip_rate_ );
		sugar_samplers_[2]->set_flip_rate( pucker_flip_rate_ );
	}
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_Suite::set_gaussian_stdev( Real const setting ) {
	gaussian_stdev_ = setting;
	if ( is_init() ) {
		chi_samplers_[1]->set_gaussian_stdev( setting );
		chi_samplers_[2]->set_gaussian_stdev( setting );
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			bb_samplers_[i]->set_gaussian_stdev( gaussian_stdev_ );
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void  RNA_MC_Suite::clear_rotamer() {
	bb_samplers_.clear();
	chi_samplers_.clear();
	sugar_samplers_.clear();
	MC_Comb::clear_rotamer();
}
///////////////////////////////////////////////////////////////////////////
} //rna
} //sampler
} //stepwise
} //protocols
