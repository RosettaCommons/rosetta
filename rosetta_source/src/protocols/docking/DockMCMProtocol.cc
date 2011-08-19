// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockMCMProtocol
/// @brief protocols that are specific to high resolution docking
/// @detailed
///		This contains the functions that create initial positions for docking
///		You can either randomize partner 1 or partner 2, spin partner 2, or
///		perform a simple perturbation.
/// 	Also contains docking mcm protocol
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
#include <protocols/moves/PackRotamersMover.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS



#include <core/pose/Pose.hh>    //JQX: dumping out pdb for testing

#include <protocols/moves/TrialMover.fwd.hh>  //JQX: add this
#include <protocols/moves/TrialMover.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.DockMCMProtocol");

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
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreFunctionCOP scorefxn_pack
) : DockingHighRes(rb_jump, scorefxn, scorefxn_pack)
{
	init();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMProtocol::DockMCMProtocol(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn
) : DockingHighRes(movable_jumps[1], scorefxn)
{
	init();
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockMCMProtocol::DockMCMProtocol(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreFunctionCOP scorefxn_pack
) : DockingHighRes(movable_jumps, scorefxn, scorefxn_pack)
{
	init();
}

//destructor
DockMCMProtocol::~DockMCMProtocol() {}

void DockMCMProtocol::set_filter( DockingHighResFilterOP filter ) { filter_ = filter; }


//clone
protocols::moves::MoverOP DockMCMProtocol::clone() const {
	return new DockMCMProtocol(*this);
}

void
DockMCMProtocol::init()
{
	moves::Mover::type( "DockMCMProtocol" );
	filter_ = new DockingHighResFilter();
}

	
	
	
	
	
	
	
	
	
	
	
/*  JQX: rewrite this apply function, try to use all the really OO concept, especial to make sure that the same mc_ object is always being used	
	
	
////////////////////////////////////////////////////////////////////////////////
/// @begin docking high resolution apply function
/// @brief
/// @detailed
///		decides what to call according to options
void DockMCMProtocol::apply( core::pose::Pose & pose )
{
	using namespace scoring;

	TR << "in DockMCMProtocol.apply" << std::endl;

	tf2()->create_and_attach_task_factory( this, pose );

	moves::PackRotamersMoverOP initial_pack = new moves::PackRotamersMover();
	initial_pack->score_function( scorefxn_pack() );
	initial_pack->task_factory( task_factory() );

	/// minimize_trial defaults to a min_tolerance of 1.0
	DockMinMoverOP minimize_trial = new DockMinMover( movable_jumps(), scorefxn() );

	dock_mcm_ = new DockMCMCycle( movable_jumps(), scorefxn(), scorefxn_pack() );
	dock_mcm_->set_task_factory( task_factory() );

	initial_pack->apply( pose );

	minimize_trial->apply( pose );

//	dock_mcm_->apply( pose );            //JQX comment this out, it seems not necessary compared to the Legacy
//	dock_mcm_->get_mc()->show_scores();  //JQX comment this out, it seems not necessary compared to the Legacy

	filter_->set_score_margin( 10.0 ); 
	if ( filter_->apply( pose ) )       
	{									
		for ( Size i=1; i<=4; ++i ) { dock_mcm_->apply( pose );
		dock_mcm_->get_mc()->show_scores();}

		pose.dump_pdb("after_4_cycles.pdb");   //JQX:  want to check the pose after 4 cycles
		exit(-1) ;  //JQX: for testing

		filter_->set_score_margin( 5.0 );  
		if ( filter_->apply( pose ) )     
		{
			for ( Size i=1; i<=45; ++i ) { dock_mcm_->apply( pose );
			dock_mcm_->get_mc()->show_scores();}
		}    
	}    
	filter_->set_score_margin( 0.0 );  
	
	// add minimize to strict tolerance 0.02
	minimize_trial->set_min_tolerance( 0.02 ); 
	minimize_trial->apply( pose );
}
	
	
*/
	
	
	
	
	
	
	
// JQX: rewrite the apply function
	
////////////////////////////////////////////////////////////////////////////////
/// @begin docking high resolution apply function
/// @brief
/// @detailed
///		decides what to call according to options
void DockMCMProtocol::apply( core::pose::Pose & pose )
{

	using namespace scoring;

	TR << "in DockMCMProtocol.apply" << std::endl;

	dock_mcm_ = new DockMCMCycle( movable_jumps(), scorefxn(), scorefxn_pack() );
	//check if we are ignoring the default docking task
	if( !ignore_default_task() ){
	 	TR << "Using the DockingTaskFactory." << std::endl;
	 	tf2()->create_and_attach_task_factory( this, pose );
	 }
	else {
	 	TR <<"The default DockingTaskFactory is being ignored." << std::endl;
	 }
	//Need a check to make sure there is a task HERE!!!
	if( task_factory() )
		dock_mcm_->set_task_factory( task_factory() );
	if( ignore_default_task() && !task_factory() ){
		TR<< "No TaskFactory given to docking, using default to design all residues." << std::endl;
	}
	//debugging
	//std::cout << "DockMCMProtocol task set: \n" <<  *(task_factory()->create_task_and_apply_taskoperations(pose)) << std::endl;

	//JQX: define the initial_pack, and this initial_pack was defined as a trial mover
	moves::PackRotamersMoverOP initial_pack = new moves::PackRotamersMover();
	initial_pack->score_function( scorefxn_pack() );
	initial_pack->task_factory( task_factory() );
	if ( dock_mcm_->get_mc()->last_accepted_pose().empty() ) { dock_mcm_->init_mc(pose); } //JQX: use the dock_mcm_'s "mc_" object
	moves::TrialMoverOP initial_pack_trial = new moves::TrialMover(initial_pack, dock_mcm_->get_mc() );
	initial_pack_trial->apply( pose );

	/// minimize_trial defaults to a min_tolerance of 0.01 in DockMinMover
	DockMinMoverOP minimize_trial = new DockMinMover( movable_jumps(), scorefxn(), dock_mcm_->get_mc() );
	minimize_trial->apply( pose );


//	filter_->set_score_margin( 10.0 );
//	if ( filter_->apply( pose ) )
//	{
		for ( Size i=1; i<=4; ++i ) {
			dock_mcm_->apply( pose );
//			dock_mcm_->get_mc()->show_scores();
		}

//		pose.dump_pdb("after_4_cycles.pdb");


//		filter_->set_score_margin( 5.0 );

		dock_mcm_->reset_cycle_index(); //JQX: reset the CycleMover index to 0


//		if ( filter_->apply( pose ) )
//		{
			for ( Size i=1; i<=45; ++i ) {
				dock_mcm_->apply( pose );
//				dock_mcm_->get_mc()->show_scores();
			}
//		}
//	}

//	filter_->set_score_margin( 0.0 );  
	

	minimize_trial->set_min_tolerance( 0.01 ); 
	minimize_trial->apply( pose );
	
	
	
	
//JQX:	get the smallest energy pose memorized in the "mc_" object
	dock_mcm_-> get_mc()->recover_low( pose );
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	

std::string
DockMCMProtocol::get_name() const {
	return "DockMCMProtocol";
}


} // namespace docking
} // namespace protocols
