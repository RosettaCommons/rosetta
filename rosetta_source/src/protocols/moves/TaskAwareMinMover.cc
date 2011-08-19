// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/TaskAwareMinMover.cc
/// @brief TaskAwareMinMover methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/moves/TaskAwareMinMover.hh>
#include <protocols/moves/TaskAwareMinMoverCreator.hh>

// Package Headers

// Project Headers
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <protocols/moves/MinMover.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;


static basic::Tracer TR( "protocols.moves.TaskAwareMinMover" );

namespace protocols {
namespace moves {

std::string
TaskAwareMinMoverCreator::keyname() const
{
	return TaskAwareMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
TaskAwareMinMoverCreator::create_mover() const {
	return new TaskAwareMinMover;
}

std::string
TaskAwareMinMoverCreator::mover_name()
{
	return "TaskAwareMinMover";
}

TaskAwareMinMover::TaskAwareMinMover()
	: Mover("TaskAwareMinMover"),
		minmover_(0),
		factory_(0),
		chi_(true),
		bb_(false)
{}

///@brief constructor with TaskFactory
TaskAwareMinMover::TaskAwareMinMover(
	protocols::moves::MinMoverOP minmover_in,
	core::pack::task::TaskFactoryCOP factory_in
) : Mover("TaskAwareMinMover"),
		minmover_(minmover_in),
		factory_(factory_in),
		chi_(true),
		bb_(false)
{
	Mover::type( "TaskAwareMinMover" );
}

TaskAwareMinMover::~TaskAwareMinMover(){}

///@details apply will extract the movemap from your minmover, modify it to include sidechain DOFs that are packable according to some TaskFactory, run the minmover with this movemap, and revert the minmover to its original movemap.
void TaskAwareMinMover::apply( core::pose::Pose & pose ){
	runtime_assert( minmover_ );
	runtime_assert( factory_ );

  using core::kinematics::MoveMapOP;
  using core::kinematics::MoveMap;

	// non-initialization failsafe
	if ( ! minmover_->movemap() ) minmover_->movemap( new MoveMap );

  //clone the MinMover's MoveMap
  MoveMapOP mm = new MoveMap( *( minmover_->movemap() ) );
  MoveMapOP mm_copy = new MoveMap( *( minmover_->movemap() ) );

  //generate task
  using core::pack::task::PackerTaskOP;
  PackerTaskOP task( factory_->create_task_and_apply_taskoperations( pose ) );

  //modify movemap by task
//  core::kinematics::modify_movemap_from_packertask( *mm, *task );
	Size const nres( task->total_residue() );

	for ( Size i(1); i <= nres; ++i ) {
		if ( task->pack_residue( i ) ) {
			// the MoveMap initializes to false for all degrees of freedom
			// this class only turns on minimization at packable dofs, it does not turn them off
			if ( chi_ ) mm->set_chi( i, chi_ );
			if ( bb_ ) mm->set_bb( i, bb_ );
		}
	}

  //pass the modified map into the MinMover
  minmover_->movemap( mm );

	//now run MinMover
	minmover_->apply( pose );

	//restore MinMover's original movemap to prevent accumulation of state
  minmover_->movemap( mm_copy );


}//apply

std::string
TaskAwareMinMover::get_name() const {
	return TaskAwareMinMoverCreator::mover_name();
}

MoverOP TaskAwareMinMover::fresh_instance() const { return new TaskAwareMinMover; }
MoverOP TaskAwareMinMover::clone() const { return new TaskAwareMinMover( *this ); }

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
TaskAwareMinMover::parse_my_tag(
	TagPtr const tag,
	DataMap & datamap,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->getName() != "TaskAwareMinMover" ) {
		TR << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("chi") ) chi_ = tag->getOption<bool>("chi");
	if ( tag->hasOption("bb") ) bb_ = tag->getOption<bool>("bb");
	minmover_ = new MinMover;
	// call to MinMover::parse_my_tag avoided here to prevent collision of chi and bb tag options
	minmover_->parse_opts( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
TaskAwareMinMover::parse_task_operations(
	TagPtr const tag,
	DataMap const & datamap,
	Filters_map const &,
	Movers_map const &,
	Pose const &
)
{
	using namespace core::pack::task;
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	factory_ = new_task_factory;
}

}//moves
}//protocols

