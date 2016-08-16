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
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>

// Package headers
// #include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/Loops.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>

// Basic headers
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_mover.refine.LoopRefineInnerCycle" );
using namespace core;

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
LoopRefineInnerCycle::LoopRefineInnerCycle() : Mover()
{
	init();
}

/// @brief copy constructor
LoopRefineInnerCycle::LoopRefineInnerCycle( LoopRefineInnerCycle const & rhs ) : Mover(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
LoopRefineInnerCycle & LoopRefineInnerCycle::operator=( LoopRefineInnerCycle const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	Mover::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
LoopRefineInnerCycle::~LoopRefineInnerCycle() {}

/// @brief Each derived class must specify its name.
std::string LoopRefineInnerCycle::get_name() const
{
	return type();
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool LoopRefineInnerCycle::reinitialize_for_new_input() const
{
	return true;
}

void LoopRefineInnerCycle::register_options()
{
	///  PUT THE LIST OF OPTIONS THAT ARE USED HERE  ///

	///  RECURSIVELY CALL REGISTER OPTIONS ON ALL MOVERS THAT THIS CLASS HAS AN OWNING_PTR TO  ///
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

// constructor with arguments
LoopRefineInnerCycle::LoopRefineInnerCycle(
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) : Mover()
{
	init( loop_mover, mc, scorefxn, tf );
}
void LoopRefineInnerCycle::setup_objects( Pose const & /* pose */ )
{
	TR << "Setting up data for " + get_name() + "." << std::endl;

	/// Perform some sanity checks to ensure the data integrity before moving forward
	using utility::excn::EXCN_Msg_Exception;

	if ( !scorefxn_ ) {
		throw EXCN_Msg_Exception( "No ScoreFunction available in " + get_name() + "." );
	}

	if ( !tf_ ) {
		throw EXCN_Msg_Exception( "No TaskFactory available in " + get_name() + "." );
	}

	if ( !mc_ ) {
		throw EXCN_Msg_Exception( "No MonteCarlo instance available in " + get_name() + "." );
	}

	if ( loop_mover_that_owns_me_.expired() ) {
		throw EXCN_Msg_Exception( "No parent LoopMover available in " + get_name() + ". This is needed to provide information on the progress of the simulation." );
	}
}

void LoopRefineInnerCycle::init()
{
	init( LoopMover_Refine_CCDAP(), NULL, NULL, NULL );
}

void LoopRefineInnerCycle::init(
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) {
	loop_mover_that_owns_me_ = loop_mover;
	mc_ = mc;
	scorefxn_ = scorefxn;
	tf_ = tf;

	type( "LoopRefineInnerCycle" );
	init_options();
}

void LoopRefineInnerCycle::init_for_equal_operator_and_copy_constructor(
	LoopRefineInnerCycle & lhs,
	LoopRefineInnerCycle const & rhs
)
{
	// copy all data members from rhs to lhs
	lhs.debug_ = rhs.debug_;
	lhs.loop_mover_that_owns_me_ = rhs.loop_mover_that_owns_me_;
	lhs.mc_ = rhs.mc_;
	lhs.scorefxn_ = rhs.scorefxn_;
	lhs.tf_ = rhs.tf_;
	lhs.movemap_ = rhs.movemap_;
}

void LoopRefineInnerCycle::init_options()
{
	using namespace basic::options;

	// Set options here.
	set_debug( option[ OptionKeys::loops::debug ].user() );
}

bool LoopRefineInnerCycle::debug() const
{
	return debug_;
}

void LoopRefineInnerCycle::set_debug( bool debug )
{
	debug_ = debug;
}

moves::MonteCarloOP LoopRefineInnerCycle::mc() const
{
	return mc_;
}

void LoopRefineInnerCycle::set_mc( moves::MonteCarloOP mc)
{
	mc_ = mc;
}

core::scoring::ScoreFunctionOP LoopRefineInnerCycle::scorefxn() const
{
	return scorefxn_;
}

void LoopRefineInnerCycle::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn)
{
	scorefxn_ = scorefxn;
}

core::pack::task::TaskFactoryOP LoopRefineInnerCycle::task_factory() const
{
	return tf_;
}

void LoopRefineInnerCycle::set_task_factory( core::pack::task::TaskFactoryOP tf )
{
	tf_ = tf;
}

core::kinematics::MoveMapOP LoopRefineInnerCycle::movemap() const
{
	// Lazily instantiate a movemap
	if ( ! movemap_ ) { movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap ); }
	return movemap_;
}

void LoopRefineInnerCycle::set_movemap( core::kinematics::MoveMapOP movemap )
{
	movemap_ = movemap;
}

Loops LoopRefineInnerCycle::get_one_random_loop() const
{
	LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
	Loops::const_iterator it( loop_mover_op->loops()->one_random_loop() );
	Loops one_loop;
	one_loop.add_loop( it );
	return one_loop;
}

LoopMover_Refine_CCDAP LoopRefineInnerCycle::loop_mover() const
{
	return loop_mover_that_owns_me_;
}

void LoopRefineInnerCycle::set_loop_mover( LoopMover_Refine_CCDAP new_owner_in_town )
{
	loop_mover_that_owns_me_ = new_owner_in_town;
}

void
LoopRefineInnerCycle::show( std::ostream & out ) const
{
	out << *this;
}

std::ostream & operator<<(std::ostream& out, LoopRefineInnerCycle const & loop_refine_inner_cycle )
{
	out << loop_refine_inner_cycle.get_name() << " is an abstract class.  Only subclasses can be used." << std::endl;
	return out;
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
