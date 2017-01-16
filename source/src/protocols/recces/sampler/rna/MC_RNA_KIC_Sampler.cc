// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.cc
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>

// Package headers
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Sugar.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh>
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>
#include <core/id/NamedAtomID.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
using namespace protocols::recces::sampler;
using namespace protocols::stepwise::sampler::screener;
using namespace protocols::stepwise::monte_carlo::mover;
using namespace protocols::stepwise::sampler::rna;

static THREAD_LOCAL basic::Tracer TR( "protocols.recces.sampler.rna.MC_RNA_KIC_Sampler" );

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

MC_RNA_KIC_Sampler::MC_RNA_KIC_Sampler(
	core::pose::PoseCOP mc_pose,
	core::Size const moving_suite,
	core::Size const chainbreak_suite,
	bool const change_ft
):
	MC_Sampler(),
	moving_suite_( moving_suite ),
	chainbreak_suite_( chainbreak_suite ),
	max_tries_( 100 ),
	verbose_( false ),
	did_close_( false ),
	cutpoint_handler_(TransientCutpointHandlerOP ( new TransientCutpointHandler( moving_suite, chainbreak_suite, change_ft )))
{
	update_pose_ = mc_pose;
}

//////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::init() {
	using namespace core::id;

	if ( is_init() ) return; // do not init more than once!

	/////// Backbone MC Sampler ////////

	TorsionID epsilon_ID( TorsionID( moving_suite_, BB, EPSILON ) );
	TorsionIDs.push_back( epsilon_ID );
	bb_samplers_.emplace_back( new MC_OneTorsion( epsilon_ID, update_pose_->torsion( epsilon_ID) ) );

	TorsionID zeta_ID( TorsionID( moving_suite_, BB, ZETA));
	TorsionIDs.push_back( zeta_ID );
	bb_samplers_.emplace_back( new MC_OneTorsion( zeta_ID, update_pose_->torsion( zeta_ID ) ) );

	TorsionID alpha1_ID( TorsionID( moving_suite_ + 1, BB, ALPHA ) );
	TorsionIDs.push_back( alpha1_ID );
	bb_samplers_.emplace_back( new MC_OneTorsion( alpha1_ID, update_pose_->torsion( alpha1_ID ) ) );

	TorsionID alpha2_ID( TorsionID( chainbreak_suite_ + 1, BB, ALPHA ) );
	TorsionIDs.push_back( alpha2_ID );
	bb_samplers_.emplace_back( new MC_OneTorsion( alpha2_ID, update_pose_->torsion( alpha2_ID ) ) );

	//Add all the pivot torsions to the list of the torsion IDs so that we have a
	//full list of the torsions that could be modified by this sampler
	TorsionIDs.emplace_back( moving_suite_ + 1, BB, BETA );
	TorsionIDs.emplace_back( moving_suite_ + 1, BB, GAMMA );
	TorsionIDs.emplace_back( moving_suite_ + 1, BB, EPSILON ); // wait a minute: chainbreak_suite_? rhiju2016
	TorsionIDs.emplace_back( moving_suite_ + 1, BB, ZETA );    // wait a minute: chainbreak_suite_? rhiju2016
	TorsionIDs.emplace_back( chainbreak_suite_ + 1, BB, BETA );
	TorsionIDs.emplace_back( chainbreak_suite_ + 1, BB, GAMMA );

	// Now find the initial torsions
	core::Real torsion;
	for ( Size i = 1; i <= TorsionIDs.size(); ++i ) {
		torsion = update_pose_->torsion(TorsionIDs[i]);
		initial_torsions_.push_back( torsion );
		angle_min_.push_back(-180);
		angle_max_.push_back(180);
	}

	//Let's not worry about the sugar rotamers for now (1/15/15)

	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->init();
	}

	////////// Make a stored loop closer /////////
	stored_loop_closer_ = RNA_KinematicCloserOP( new RNA_KinematicCloser(	*update_pose_, moving_suite_, chainbreak_suite_ ) );
	stored_loop_closer_->set_calculate_jacobians( true );

	////////// Loop Closer //////////
	loop_closer_ = RNA_KinematicCloserOP( new RNA_KinematicCloser( *update_pose_, moving_suite_, chainbreak_suite_ ) );
	loop_closer_->set_calculate_jacobians( true );
	loop_closer_->set_verbose( verbose_ );

	set_init( true );
}

////////////////////////////////////////////////////////////////////////
void
MC_RNA_KIC_Sampler::get_next_solutions( pose::Pose const & pose )
{
	Pose pose_copy = pose;
	//Need to update the stored angles in the samplers in case they were changed in a different sampler
	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->set_angle( pose_copy.torsion( TorsionIDs[i] ) ) ;
	}

	//First need to make a cutpoint
	//Then find the solutions for the closed pose
	cutpoint_handler_->put_in_cutpoints( pose_copy );
	stored_loop_closer_->init(pose_copy, pose);
	stored_jacobians_ = stored_loop_closer_->get_all_jacobians();

	// Then, find any solutions for pose with perturbed 'driver' torsions.
	for ( Size i = 1; i <= max_tries_; ++i ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			++( *bb_samplers_[i]);
			bb_samplers_[i]->apply( pose_copy );
		}

		loop_closer_->init(pose_copy, pose);
		//if there are 0 solutions, try again
		if ( !loop_closer_->not_end() ) {
			continue;
		}
		cutpoint_handler_->take_out_cutpoints( pose_copy );
		did_close_ =  true ;
		return;
	}

	cutpoint_handler_->take_out_cutpoints( pose_copy );
	did_close_ = false ;

	TR.Debug << "Chain not closable after " << max_tries_ << " tries!" << std::endl;

}
////////////////////////////////////////////////////////////////////////
void
MC_RNA_KIC_Sampler::choose_solution() {
	if ( did_close_ ) { // loop closed after bb perturbations in get_next_solutions()
		current_jacobians_ = loop_closer_->get_all_jacobians();
		all_jacobians_ = current_jacobians_; // perturbed solutions
		all_jacobians_.append( stored_jacobians_ ); // closing solutions on stored pose.
		jacobian_sampler_ = numeric::random::WeightedSampler( all_jacobians_ );
		solution_ = jacobian_sampler_.random_sample( numeric::random::rg().uniform() );
	} else {
		all_jacobians_ = stored_jacobians_;
		jacobian_sampler_ = numeric::random::WeightedSampler( all_jacobians_ );
		solution_ = jacobian_sampler_.random_sample( numeric::random::rg().uniform() );
	}
}
////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::operator++() {
	next( *update_pose_ );
}
////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::next( Pose const & pose ) {
	get_next_solutions( pose );
	choose_solution();
}
///////////////////////////////////////////////////////////////////////////
bool MC_RNA_KIC_Sampler::check_moved() const {
	return found_move_;
}
///////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::apply( pose::Pose & pose_in ) {
	runtime_assert( is_init() );
	Pose pose = pose_in;
	cutpoint_handler_->put_in_cutpoints( pose );
	used_current_solution_ = false;
	found_move_ = false;

	//// apply solution_ here!
	if ( did_close_ ) {
		if ( solution_ <= current_jacobians_.size() ) {
			for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
				bb_samplers_[i]->apply( pose );
			}
			loop_closer_->apply( pose, solution_ );
			used_current_solution_ = true;
			found_move_ = true;
		} else {
			Real calculated_jacobian = get_jacobian( pose );
			stored_loop_closer_->apply( pose, (solution_ - current_jacobians_.size()) );
			//Check whether this is actually a new pose
			//This isn't really a good way to check if the pose is the same...change this!!
			Real picked_jacobian = stored_jacobians_[ solution_ - current_jacobians_.size()];
			if ( std::abs(calculated_jacobian - picked_jacobian) > 0.0000000001 ) { found_move_ = true; }
			//used_current solution remains false, so don't update stored angles, don't update stored loop closer
		}
	} else {
		Real calculated_jacobian = get_jacobian( pose );
		stored_loop_closer_->apply( pose, solution_ );
		//Check whether this is actually a new pose
		//This isn't really a good way to check if the pose is the same...change this!!
		Real picked_jacobian = stored_jacobians_[ solution_ - current_jacobians_.size()];
		if ( std::abs(calculated_jacobian - picked_jacobian) > 0.0000000001 ) { found_move_ = true; }
		//used_current solution remains false, so don't update stored angles, don't update stored loop closer
	}

	cutpoint_handler_->take_out_cutpoints( pose );
	if ( found_move_ ) {
		if ( !( check_angles_in_range( pose ) ) ) {
			found_move_ = false;
			return;
		}
	}
	pose_in = pose;
}
///////////////////////////////////////////////////////////////////////////
Real MC_RNA_KIC_Sampler::get_jacobian( pose::Pose & pose ) {
	return loop_closer_->get_jacobian(pose);
}
///////////////////////////////////////////////////////////////////////////
Real MC_RNA_KIC_Sampler::vector_sum( utility::vector1< core::Real > const & vector ) {
	sum_ = 0;
	for ( Size i = 1; i<=vector.size(); ++i ) {
		sum_ += vector[i];
	}
	return sum_;
}
///////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::update() {
	runtime_assert( is_init() );
	if ( used_current_solution_ ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			bb_samplers_[i]->update();
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::set_gaussian_stdev( Real const setting ) {
	gaussian_stdev_ = setting;
	if ( is_init() ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			bb_samplers_[i]->set_gaussian_stdev( gaussian_stdev_ );
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void MC_RNA_KIC_Sampler::set_angle_range_from_init_torsions( core::Real const range ) {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->set_angle_range( initial_torsions_[i] - range, initial_torsions_[i] + range );
	}
	for ( Size i = bb_samplers_.size() + 1; i <= initial_torsions_.size(); ++i ) {
		angle_min_[i] = initial_torsions_[i] - range;
		angle_max_[i] = initial_torsions_[i] + range;
	}
}
/////////////////////////////////////////////////////////////////////////
bool MC_RNA_KIC_Sampler::check_angles_in_range( const Pose & pose ) {
	// Just need to check the torsions for the pivot angles, the ones for the driver angles
	// are already guaranteed to be in the correct range by the set_angle_range fn
	// (see set_angle_range_from_init_torsions)
	for ( Size i = bb_samplers_.size() + 1; i <= TorsionIDs.size(); i++ ) {
		if ( angle_max_[i] - angle_min_[i] >= 360 ) continue;
		Real angle( pose.torsion( TorsionIDs[i] ) );
		while ( angle < angle_min_[i] ) angle += 360;
		while ( angle >= angle_max_[i] ) angle -= 360;
		if ( angle >= angle_min_[i] ) continue;
		return false;
	}
	return true;
}
///////////////////////////////////////////////////////////////////////////
void
MC_RNA_KIC_Sampler::show( std::ostream & out, Size const indent ) const {
	for ( Size n = 1; n <= indent; n++ ) out << ' ';
	out << ( "MC_RNA_KIC_Sampler moving_suite:" + utility::to_string(moving_suite_) + " chainbreak_suite:" +
		utility::to_string(chainbreak_suite_) );
}
/////////////////////////////////////////////////////////////////////////
/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
MC_SamplerOP
MC_RNA_KIC_Sampler::find( core::id::TorsionID const & torsion_id )
{
	if ( get_suite_torsion_ids( moving_suite_ ).has_value( torsion_id ) ||
			get_suite_torsion_ids( chainbreak_suite_ ).has_value( torsion_id ) ) {
		return std::dynamic_pointer_cast< MC_Sampler >( shared_from_this() );
	}
	return 0;
}
/////////////////////////////////////////////////////////////////////////
} //rna
} //sampler
} //recces
} //protocols
