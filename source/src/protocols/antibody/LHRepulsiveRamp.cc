// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/LHRepulsiveRamp.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#include <protocols/antibody/LHRepulsiveRamp.hh>
#include <basic/Tracer.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/docking/DockMCMCycle.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.LHRepulsiveRamp" );


using namespace core;
namespace protocols {
namespace antibody {


// default constructor
LHRepulsiveRamp::LHRepulsiveRamp() : Mover() {}


LHRepulsiveRamp::LHRepulsiveRamp(  docking::DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP dock_scorefxn,
	core::scoring::ScoreFunctionCOP pack_scorefxn ) : Mover() {
	user_defined_ = true;
	jump_ = movable_jumps;
	dock_scorefxn_ = dock_scorefxn->clone();
	pack_scorefxn_ = pack_scorefxn->clone();

	init();
}


// default destructor
LHRepulsiveRamp::~LHRepulsiveRamp() {}

//clone
protocols::moves::MoverOP LHRepulsiveRamp::clone() const {
	return( protocols::moves::MoverOP( new LHRepulsiveRamp() ) );
}


void LHRepulsiveRamp::init( ) {
	set_default();


	//JQX: Jeff wants this repulsive ramping mover to be more general, therefore, no default
	//      tf, movemap, and scorefxn are set up here, to avoid any unnecessary confusion.
}


void LHRepulsiveRamp::set_default() {
	benchmark_       = false;
	sc_min_ =false;
	rt_min_ =false;


	rep_ramp_cycles_ = 3 ;
	rot_mag_         = 2.0 ;
	trans_mag_       = 0.1 ;
	num_repeats_     = 4;

	if ( !user_defined_ ) {}

}


///////////////////////////////////////////////////////////////////////////
///
/// @brief ramping up the fullatom repulsive weight slowly to allow the
///        partners to relieve clashes and make way for each other
///
/// @details This routine is specially targetted to the coupled
///           optimization of docking partners and the loop region.  The
///           loop modelling & all previous  steps  involve mainly
///           centroid  mode .On switching  on fullatom mode, one is bound
///           to end up with clashes.To relieve the clashes, it is
///           essential to slowly  dial up the  repulsive weight instead of
///           turning it on to the maximum  value in one single step
///
/// @param[in] input pose which is assumed to have a docking fold tree
///
/// @global_read fa_rep : fullatom repulsive weight
///
/// @remarks A particular portion is  commented out,which can be
///          uncommented if one  uses a  low  resolution  homology  model.
///          Check details in the beginning of the commented out region
///
///////////////////////////////////////////////////////////////////////////

void LHRepulsiveRamp::apply( pose::Pose & pose ) {
	TR<<"start apply function ..."<<std::endl;

	//JQX: the fold_tree and variants should have been set up in "pose"
	//     executed outside of this class


	// add scores to map
	( *dock_scorefxn_ )( pose );

	// dampen fa_rep weight
	core::Real rep_weight_max = dock_scorefxn_->get_weight( core::scoring::fa_rep );

	if ( benchmark_ ) {
		rep_ramp_cycles_ = 1;
		num_repeats_ = 1;
	}


	core::Real rep_ramp_step = (rep_weight_max - 0.02) / core::Real(rep_ramp_cycles_-1);
	core::scoring::ScoreFunctionOP temp_dock_scorefxn = dock_scorefxn_->clone();

	for ( Size i = 1; i <= rep_ramp_cycles_; i++ ) {
		core::Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
		TR<<"   repulsive ramp cycle "<<i<<":     rep_weight = "<<rep_weight<<std::endl;
		temp_dock_scorefxn->set_weight( core::scoring::fa_rep, rep_weight );


		docking::DockMCMCycleOP dockmcm_cyclemover( new docking::DockMCMCycle( jump_, temp_dock_scorefxn, pack_scorefxn_ ) );
		//TODO: print scoring function in apply and move "new" out
		dockmcm_cyclemover->set_rot_magnitude(rot_mag_);
		dockmcm_cyclemover->set_task_factory(tf_);
		dockmcm_cyclemover->set_move_map(movemap_);
		if ( sc_min_ ) {
			dockmcm_cyclemover->set_scmin(true);
		}
		if ( sc_min_ ) {
			dockmcm_cyclemover->set_scmin(true);
		}

		for ( Size j=1; j<=num_repeats_; j++ ) {
			dockmcm_cyclemover -> apply(pose);
			TR<<"       doing rb_mover_min_trial in the DockMCMCycle  ...   "<<j<<std::endl;
		}
		dockmcm_cyclemover -> reset_cycle_index(); //JQX: only do the rb_mover_min_trial (index<7)
		dockmcm_cyclemover -> get_mc()->recover_low( pose ); //chose the pose having the lowest score

	}

	TR<<"finish apply function !!!"<<std::endl;

}


std::string LHRepulsiveRamp::get_name() const {
	return "LHRepulsiveRamp";
}

void LHRepulsiveRamp::set_task_factory(pack::task::TaskFactoryCOP tf) {
	tf_ = pack::task::TaskFactoryOP( new pack::task::TaskFactory(*tf) );
}

void LHRepulsiveRamp::set_move_map(kinematics::MoveMapCOP movemap) {
	movemap_ = kinematics::MoveMapOP( new kinematics::MoveMap(*movemap) );
}


void LHRepulsiveRamp::set_dock_jump(docking::DockJumps jump) {
	jump_ = jump;
}


} // namespace antibody
} // namespace protocols

