// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockMCMProtocol
/// @brief protocols that are specific to high resolution docking
/// @details
///  This contains the functions that create initial positions for docking
///  You can either randomize partner 1 or partner 2, spin partner 2, or
///  perform a simple perturbation.
///  Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn

// Unit Headers
#include <protocols/docking/DockMCMProtocol.hh>

// Package Headers
#include <protocols/docking/DockFilters.hh>
#include <protocols/docking/DockMinMover.hh>
#include <protocols/docking/DockMCMCycle.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/kinematics/FoldTree.hh>
// Utility Headers
#include <basic/Tracer.hh>

#include <protocols/jd2/util.hh>

#include <core/pose/Pose.hh>    //JQX: dumping out pdb for testing

#include <protocols/moves/TrialMover.fwd.hh>  //JQX: add this
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <ObjexxFCL/format.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <protocols/docking/DockTaskFactory.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockMCMProtocol", basic::t_info );

//     originally from dock_structure.cc Jeff Gray April 2001

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockMCMProtocol::DockMCMProtocol() : DockingHighRes( 1 )
{
	init();
}

// constructor with arguments
// only one movable jump
DockMCMProtocol::DockMCMProtocol(
	core::Size const rb_jump
) : DockingHighRes(rb_jump)
{
	init();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMProtocol::DockMCMProtocol(
	core::Size const rb_jump,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : DockingHighRes(rb_jump, scorefxn, scorefxn_pack)
{
	init();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMProtocol::DockMCMProtocol(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionOP scorefxn
) : DockingHighRes(movable_jumps[1], scorefxn)
{
	init();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMProtocol::DockMCMProtocol(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : DockingHighRes(movable_jumps, scorefxn, scorefxn_pack)
{
	init();
}

//destructor
DockMCMProtocol::~DockMCMProtocol() {}

void DockMCMProtocol::set_filter( DockingHighResFilterOP filter ) { filter_ = filter; }


//clone
protocols::moves::MoverOP DockMCMProtocol::clone() const {
	return protocols::moves::MoverOP( new DockMCMProtocol(*this) );
}

void
DockMCMProtocol::init()
{
	moves::Mover::type( "DockMCMProtocol" );
	filter_ = DockingHighResFilterOP( new DockingHighResFilter() );
	movemap_reset_ = false;
	num_of_first_cycle_=4;
	num_of_second_cycle_=45;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ OptionKeys::docking::dock_mcm_first_cycles ].user() ) {
		set_first_cycle(option[ OptionKeys::docking::dock_mcm_first_cycles ]() );
	}
	if ( option[ OptionKeys::docking::dock_mcm_second_cycles ].user() ) {
		set_second_cycle(option[ OptionKeys::docking::dock_mcm_second_cycles ]() );
	}


}

// JQX: rewrite the apply function

////////////////////////////////////////////////////////////////////////////////
/// @brief
/// @details
///  decides what to call according to options
void DockMCMProtocol::apply( core::pose::Pose& pose )
{

	using namespace scoring;

	TR.Info << "in DockMCMProtocol.apply" << std::endl;
	TR.Debug << "fold-tree in DockMCM: " << pose.fold_tree() << std::endl;
	(*scorefxn())( pose );
	jd2::write_score_tracer( pose, "DockMCM_start" );

	dock_mcm_ = DockMCMCycleOP( new DockMCMCycle( movable_jumps(), scorefxn(), scorefxn_pack() ) );
	if ( movemap_reset_ ) {
		dock_mcm_->set_move_map(movemap_);
	}
	//check if we are ignoring the default docking task
	if ( !ignore_default_task() ) {
		TR << "Using the DockingTaskFactory." << std::endl;
		tf2()->create_and_attach_task_factory( this, pose );
	} else {
		TR <<"The default DockingTaskFactory is being ignored." << std::endl;
	}
	//Need a check to make sure there is a task HERE!!!
	if ( task_factory() != 0 ) {
		dock_mcm_->set_task_factory( task_factory() );
	}
	//exit if you chose to ignore the default task but didn't provide one of your own.
	if ( ignore_default_task() && !task_factory() ) {
		utility_exit_with_message( "Exiting DockMCMProtocol you chose to ignore_default_task but no alternate task was given to docking, " );
	}

	//JQX: define the initial_pack, and this initial_pack was defined as a trial mover
	protocols::simple_moves::PackRotamersMoverOP initial_pack( new protocols::simple_moves::PackRotamersMover() );
	initial_pack->score_function( scorefxn_pack() );
	initial_pack->task_factory( task_factory() );
	if ( dock_mcm_->get_mc()->last_accepted_pose().empty() ) { dock_mcm_->init_mc(pose); } //JQX: use the dock_mcm_'s "mc_" object
	moves::TrialMoverOP initial_pack_trial( new moves::TrialMover(initial_pack, dock_mcm_->get_mc() ) );



	//JQX: rt_min and sc_min options were ignored in the extreme coding week, put them back now
	protocols::moves::SequenceMoverOP initial_repack_sequence( new protocols::moves::SequenceMover() );
	initial_repack_sequence->add_mover(initial_pack_trial);

	if ( rt_min() ) {
		simple_moves::RotamerTrialsMinMoverOP rtmin( new simple_moves::RotamerTrialsMinMover( scorefxn_pack(), task_factory() ) );
		moves::TrialMoverOP rtmin_trial( new moves::TrialMover( rtmin, dock_mcm_->get_mc() ) );
		initial_repack_sequence->add_mover(rtmin_trial);
		dock_mcm_->set_rtmin(true);
	}
	if ( sc_min() ) {
		docking::SidechainMinMoverOP scmin_mover( new docking::SidechainMinMover( scorefxn_pack(), task_factory() ) );
		moves::TrialMoverOP scmin_trial( new moves::TrialMover( scmin_mover, dock_mcm_->get_mc() ) );
		initial_repack_sequence->add_mover(scmin_trial);
		dock_mcm_->set_scmin(true);
	}


	initial_repack_sequence->apply( pose );



	jd2::write_score_tracer( pose, "DockMCM_pack_trialed" );

	/// minimize_trial defaults to a min_tolerance of 0.01 in DockMinMover
	DockMinMoverOP minimize_trial( new DockMinMover( movable_jumps(), scorefxn(), dock_mcm_->get_mc() ) );
	minimize_trial->apply( pose );
	jd2::write_score_tracer( pose, "DockMCM_minimized" );

	// filter_->set_score_margin( 10.0 );
	// if ( filter_->apply( pose ) )
	// {
	for ( Size i=1; i<=num_of_first_cycle_; ++i ) {
		dock_mcm_->apply( pose );
		jd2::write_score_tracer( pose, "DockMCM_cycle_"+ObjexxFCL::format::I(1, i ) );
	}

	//  filter_->set_score_margin( 5.0 );

	dock_mcm_->reset_cycle_index(); //JQX: reset the CycleMover index to 0


	//  if ( filter_->apply( pose ) )
	//  {
	for ( Size i=1; i<=num_of_second_cycle_; ++i ) {
		dock_mcm_->apply( pose );
		jd2::write_score_tracer( pose, "DockMCM_cycle_"+ObjexxFCL::format::I(1, i+4 ) );
	}
	//  }
	// }

	// filter_->set_score_margin( 0.0 );


	minimize_trial->set_min_tolerance( 0.01 );
	minimize_trial->apply( pose );

	//JQX: get the smallest energy pose memorized in the "mc_" object
	dock_mcm_-> get_mc()->recover_low( pose );

	jd2::write_score_tracer( pose, "DockMCM_final" );

}


void DockMCMProtocol::set_move_map(core::kinematics::MoveMapOP movemap){
	movemap_reset_ = true;
	movemap_ = movemap;
}

void DockMCMProtocol::set_first_cycle(Size const & num){
	num_of_first_cycle_ = num;
}
void DockMCMProtocol::set_second_cycle(Size const & num){
	num_of_second_cycle_ = num;
}

core::scoring::ScoreFunctionCOP DockMCMProtocol::scorefxn_docking() const {
	return scorefxn();
}

core::scoring::ScoreFunctionCOP DockMCMProtocol::scorefxn_packing() const {
	return scorefxn_pack();
}

std::string
DockMCMProtocol::get_name() const {
	return "DockMCMProtocol";
}

std::ostream & operator<<(std::ostream& out, const DockMCMProtocol & dmp )
{
	// Display the state of the filters (on or off)
	out << "High Resolution Filter: " << ( (dmp.filter_) ? "on" : "off" ) << std::endl;
	out << "Docking Scorefunction:  " << dmp.scorefxn_docking()->get_name() << std::endl;
	out << "Packing Scorefunction:  " << dmp.scorefxn_packing()->get_name() << std::endl;
	out << "Movemap: ";
	if ( dmp.movemap_ != 0 ) { out << std::endl; dmp.movemap_->show(); }
	else { out << "none"; }
	return out;
}

} // namespace docking
} // namespace protocols
