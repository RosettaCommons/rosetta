// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockingLowResEnsemble
/// @brief protocols that are specific to low resolution ensemble docking
/// @details
/// @author Daisuke Kuroda

#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingLowResEnsemble.hh>
#include <protocols/docking/DockingEnsemble.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>


#include <core/pose/Pose.hh>

#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>

// Utility Headers
#include <utility/tools/make_vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>
#include <cmath>

//Utility Headers

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.docking.DockingLowResEnsemble" );

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockingLowResEnsemble::DockingLowResEnsemble() : DockingLowRes()
{
	type( "DockingLowResEnsemble" );
}

// constructor with arguments
DockingLowResEnsemble::DockingLowResEnsemble(
	core::scoring::ScoreFunctionCOP scorefxn,
	DockJumps const movable_jumps
) : DockingLowRes( scorefxn, movable_jumps )
{
	type( "DockingLowResEnsemble" );
}

//destructor
DockingLowResEnsemble::~DockingLowResEnsemble() = default;


protocols::moves::MoverOP DockingLowResEnsemble::clone() const {
	return protocols::moves::MoverOP( new DockingLowResEnsemble(*this) );
}

core::Size DockingLowResEnsemble::get_current_conformer_ensemble1() { return current_conformer_ensemble1_; }

core::Size DockingLowResEnsemble::get_current_conformer_ensemble2() { return current_conformer_ensemble2_; }

void DockingLowResEnsemble::set_ensemble1( DockingEnsembleOP ensemble1 )
{
	ensemble1_mover_ = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble1, true ) );
}

void DockingLowResEnsemble::set_ensemble2( DockingEnsembleOP ensemble2 )
{
	ensemble2_mover_ = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble2, true ) );
}

void DockingLowResEnsemble::ensemble_defaults(){

	set_inner_cycles( 50 ); //change this
	set_outer_cycles( 10 ); //change this
	current_outer_cycle_ = 0;

	rb_trial_wt_ = 0.5;
	ens1_trial_wt_ = 0.25;
	ens2_trial_wt_ = 0.25;

	current_conformer_ensemble1_ = 1;
	current_conformer_ensemble2_ = 1;

	reset_ensemble_trial_counters();
}

void DockingLowResEnsemble::initialize_movers(){
	using namespace moves;

	get_rb_mover()->rot_magnitude( get_rot_magnitude() );
	get_rb_mover()->trans_magnitude( get_trans_magnitude() );

	set_rb_trial( TrialMoverOP ( new TrialMover( get_rb_seq(), get_mc() ) ) );

	ensemble1_trial_ = TrialMoverOP( new TrialMover( ensemble1_seq_, get_mc() ) );
	ensemble1_trial_->keep_stats_type( accept_reject );

	ensemble2_trial_ = TrialMoverOP( new TrialMover( ensemble2_seq_, get_mc() ) );
	ensemble2_trial_->keep_stats_type( accept_reject );
}

void DockingLowResEnsemble::reset_ensemble_trial_counters() {
	accept_rate_ens1_ = 0.0;
	num_trials_ens1_ = 0.0;
	num_accepted_ens1_ = 0.0;

	accept_rate_ens2_ = 0.0;
	num_trials_ens2_ = 0.0;
	num_accepted_ens2_ = 0.0;
}

void DockingLowResEnsemble::calc_ensemble_accept_rate_(){
	accept_rate_ens1_ = num_accepted_ens1_/num_trials_ens1_;
	accept_rate_ens2_ = num_accepted_ens2_/num_trials_ens2_;
}

void DockingLowResEnsemble::ensemble1_trial( core::pose::Pose & pose ){
	using namespace moves;


	int current_num_accepted = ensemble1_trial_->num_accepts();

	ensemble1_trial_->apply( pose );
	++ num_trials_ens1_;
	if ( ensemble1_trial_->num_accepts() > current_num_accepted ) {
		TR << "Accepted Conformer Switch for Ensemble 1" << std::endl;
		current_conformer_ensemble1_ = ensemble1_mover_->get_ensemble()->get_current_confnum();
		//TR << "Current Conformer Ensemble 1: " << current_conformer_ensemble1_ <<std::endl;
		++ num_accepted_ens1_;

		if ( get_view_in_pymol() ) {
			pymol_color_change( pose );
		}

	} else {
		TR << "Rejected Conformer Switch for Ensemble 1" << std::endl;
		//TR << "Current Conformer Ensemble 1: " << current_conformer_ensemble1_ <<std::endl;
	}

}

void DockingLowResEnsemble::ensemble2_trial( core::pose::Pose & pose ){
	using namespace moves;


	int current_num_accepted = ensemble2_trial_->num_accepts();

	ensemble2_trial_->apply( pose );
	++ num_trials_ens2_;
	if ( ensemble2_trial_->num_accepts() > current_num_accepted ) {
		TR << "Accepted Conformer Switch for Ensemble 2" << std::endl;
		current_conformer_ensemble2_ = ensemble2_mover_->get_ensemble()->get_current_confnum();
		//TR << "Current Conformer Ensemble 2: " << current_conformer_ensemble2_ <<std::endl;
		++ num_accepted_ens2_;

		if ( get_view_in_pymol() ) {
			pymol_color_change( pose );
		}

	} else {
		TR << "Rejected Conformer Switch for Ensemble 2" << std::endl;
		//TR << "Current Conformer Ensemble 2: " << current_conformer_ensemble2_ <<std::endl;
	}

}

void DockingLowResEnsemble::lowres_inner_cycle( core::pose::Pose & pose ) {

	adaptive_inner_cycle( pose );

	if ( current_outer_cycle_ < 8 ) {
		++ current_outer_cycle_;
	} else if ( current_outer_cycle_ == 8 ) {
		try_all_conformers( pose );
		++ current_outer_cycle_;
	} else {
		try_all_conformers( pose );
		// This is necessary to set the cacheable data in the pose to the correct conformer number 1
		ensemble1_mover_->get_ensemble()->set_current_confnum( get_current_conformer_ensemble1() );
		ensemble2_mover_->get_ensemble()->set_current_confnum( get_current_conformer_ensemble2() );
		current_outer_cycle_ = 0;
	}

	if ( accept_rate_ens1_ < 0.3 ) {
		ens1_trial_wt_ *= 1.25;
	} else {
		ens1_trial_wt_ *= 0.75;
	}

	if ( accept_rate_ens2_ < 0.3 ) {
		ens2_trial_wt_ *= 1.25;
	} else {
		ens2_trial_wt_ *= 0.75;
	}

}


void DockingLowResEnsemble::adaptive_inner_cycle( core::pose::Pose & pose ) {

	using namespace moves;

	initialize_movers();

	reset_rb_trial_counters();
	reset_ensemble_trial_counters();

	core::Real weight_sum = rb_trial_wt_ + ens1_trial_wt_ + ens2_trial_wt_;

	core::Real movechoice;

	// number of cycles in number of inner cycles specified by the user (default 50) + number of ens1 trials decided + number of ens2 trials decided
	core::Size num_inner_cycles = (core::Size)ceil(get_inner_cycles() * ( 1 + ( ens1_trial_wt_/rb_trial_wt_ ) + ( ens2_trial_wt_/rb_trial_wt_ ) ));

	for ( core::Size i = 1; i < num_inner_cycles ; ++i ) {

		movechoice = numeric::random::rg().uniform() * weight_sum;

		if ( movechoice < rb_trial_wt_ ) {
			rigid_body_trial( pose );
		} else if ( movechoice < rb_trial_wt_ + ens1_trial_wt_ ) {
			ensemble1_trial( pose );
		} else {
			ensemble2_trial( pose );
		}

	}

	pose = get_mc()->lowest_score_pose();
	calc_accept_rate();
	calc_ensemble_accept_rate_();
	TR << "Rigid body perturbation acceptance rate is: " << get_accept_rate() << std::endl;
	TR << "Conformer switch acceptance rate for ensemble 1 is: " << accept_rate_ens1_ << std::endl;
	TR << "Conformer switch acceptance rate for ensemble 2 is: " << accept_rate_ens2_ << std::endl;
	get_mc()->reset( pose );
}


void DockingLowResEnsemble::try_all_conformers( core::pose::Pose & pose ) {

	using namespace moves;

	initialize_movers();

	for ( core::Size j=1; j<=ensemble1_mover_->get_ensemble()->size(); ++j ) {
		ensemble1_mover_->set_specific_conformer( j );
		ensemble1_trial( pose );
	}

	for ( core::Size k=1; k<=ensemble2_mover_->get_ensemble()->size(); ++k ) {
		ensemble2_mover_->set_specific_conformer( k );
		ensemble2_trial( pose );
	}

	ensemble1_mover_->set_specific_conformer( 0 ); //resetting to random conformer
	ensemble2_mover_->set_specific_conformer( 0 );

	pose = get_mc()->lowest_score_pose();
	get_mc()->reset( pose );
}


void DockingLowResEnsemble::finalize_setup( core::pose::Pose & pose ){
	using namespace moves;

	DockingLowRes::finalize_setup( pose );

	ensemble_defaults(); // cannot override set_defaults() because it is called from the constructor of the base class

	ensemble1_seq_ = protocols::moves::SequenceMoverOP( new SequenceMover );
	ensemble1_seq_->add_mover( ensemble1_mover_ );

	ensemble2_seq_ = protocols::moves::SequenceMoverOP( new SequenceMover );
	ensemble2_seq_->add_mover( ensemble2_mover_ );

	if ( get_view_in_pymol() ) {
		ensemble1_seq_->add_mover( get_pymol_mover() );
		ensemble2_seq_->add_mover( get_pymol_mover() );
	}

}


void DockingLowResEnsemble::show( std::ostream & out ) const
{
	DockingLowRes::show( out );

	using namespace ObjexxFCL::format;
	std::string line_marker = "///";
	out << line_marker << " Ensemble 1: " << ( ( ensemble1_mover_ ) ? ( "on" ) : ( "off " ) );
	out << space( 59 ) << line_marker << std::endl;
	out << line_marker << " Ensemble 2: " << ( ( ensemble2_mover_ ) ? ( "on" ) : ( "off " ) );
	out << space( 59 ) << line_marker << std::endl;
}

} // namespace docking
} // namespace protocols
