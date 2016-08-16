// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/LoopRefineInnerCycle.cc
/// @brief Abstract class to define interface for all types of "inner cycle" operations used for loop refinement.
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

// Unit headers
#include <protocols/loops/loop_mover/refine/RepackTrial.hh>
#include <protocols/loops/loop_mover/refine/RepackTrialCreator.hh>

// Package headers
//#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project headers
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_mover.refine.RepackTrial" );
using namespace core;

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
RepackTrial::RepackTrial() : LoopRefineInnerCycle()
{
	init();
}

/// @brief copy constructor
RepackTrial::RepackTrial( RepackTrial const & rhs ) : LoopRefineInnerCycle(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
RepackTrial & RepackTrial::operator=( RepackTrial const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	LoopRefineInnerCycle::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
RepackTrial::~RepackTrial() {}

/// @brief Each derived class must specify its name.
std::string RepackTrial::get_name() const
{
	return type();
}

//@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
RepackTrial::clone() const
{
	return protocols::moves::MoverOP( new RepackTrial( *this ) );
}

/// @brief fresh_instance returns a default-constructed object for JD2
protocols::moves::MoverOP
RepackTrial::fresh_instance() const
{
	return protocols::moves::MoverOP( new RepackTrial() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool RepackTrial::reinitialize_for_new_input() const
{
	return true;
}

void RepackTrial::register_options()
{
	///  PUT THE LIST OF OPTIONS THAT ARE USED HERE  ///

	///  RECURSIVELY CALL REGISTER OPTIONS ON ALL MOVERS THAT THIS CLASS HAS AN OWNING_PTR TO  ///
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

// constructor with arguments
RepackTrial::RepackTrial(
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) : LoopRefineInnerCycle( loop_mover, mc, scorefxn, tf )
{
	init();
}

void RepackTrial::apply( Pose & pose )
{
	// TR << "Beginning apply function of " + get_name() + "." << std::endl;

	setup_objects( pose );

	// show( TR );

	//main_repack_trial
	LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
	if ( (loop_mover_op->current_cycle_number() % loop_mover_op->repack_period() ) == 0 ||
			loop_mover_op->current_cycle_number() == loop_mover_op->inner_cycles() ) {
		// repack trial
		pack::task::PackerTaskOP task = task_factory()->create_task_and_apply_taskoperations( pose );
		task->set_bump_check( true );

		core::pack::pack_rotamers( pose, *scorefxn(), task );
		std::string move_type = "repack";
		mc()->boltzmann( pose, move_type );
		mc()->show_scores();
	}
}


void RepackTrial::setup_objects( Pose const & pose )
{
	// TR << "Setting up data for " + get_name() + "." << std::endl;

	LoopRefineInnerCycle::setup_objects( pose );
}

void RepackTrial::init()
{
	type( "RepackTrial" );
	init_options();
}

void RepackTrial::init_for_equal_operator_and_copy_constructor(
	RepackTrial & /* lhs */,
	RepackTrial const & /* rhs */
)
{
	// copy all data members from rhs to lhs
}

void RepackTrial::init_options()
{
	/* UNCOMMENT WHEN THERE ARE ACTUALLY OPTIONS TO PROCESS
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	*/
	// Set options here.
}

void
RepackTrial::show( std::ostream & out ) const
{
	out << *this;
}

std::ostream & operator<<(std::ostream& out, RepackTrial const & repack_trial )
{
	out << repack_trial.get_name() << " is an awesome class." << std::endl;
	return out;
}

RepackTrialCreator::~RepackTrialCreator() {}

moves::MoverOP RepackTrialCreator::create_mover() const {
	return moves::MoverOP( new RepackTrial() );
}

std::string RepackTrialCreator::keyname() const {
	return "RepackTrial";
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
