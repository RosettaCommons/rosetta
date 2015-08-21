// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/TaskAwareMinMover.cc
/// @brief TaskAwareMinMover methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMoverCreator.hh>

// Package Headers

// Project Headers
#include <core/kinematics/MoveMap.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;


static thread_local basic::Tracer TR( "protocols.simple_moves.TaskAwareMinMover" );

namespace protocols {
namespace simple_moves {

std::string
TaskAwareMinMoverCreator::keyname() const
{
	return TaskAwareMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
TaskAwareMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TaskAwareMinMover );
}

std::string
TaskAwareMinMoverCreator::mover_name()
{
	return "TaskAwareMinMover";
}

TaskAwareMinMover::TaskAwareMinMover()
: protocols::moves::Mover("TaskAwareMinMover"),
	minmover_(/* 0 */),
	factory_(/* 0 */),
	chi_(true),
	bb_(false)
{}

/// @brief constructor with TaskFactory
TaskAwareMinMover::TaskAwareMinMover(
	protocols::simple_moves::MinMoverOP minmover_in,
	core::pack::task::TaskFactoryCOP factory_in
) : protocols::moves::Mover("TaskAwareMinMover"),
	minmover_(minmover_in),
	factory_(factory_in),
	chi_(true),
	bb_(false),
	jump_(false)
{
	protocols::moves::Mover::type( "TaskAwareMinMover" );
}

TaskAwareMinMover::~TaskAwareMinMover(){}

/// @details apply will extract the movemap from your minmover, modify it to include sidechain DOFs that are packable according to some TaskFactory, run the minmover with this movemap, and revert the minmover to its original movemap.
void TaskAwareMinMover::apply( core::pose::Pose & pose ){
	runtime_assert( minmover_ != 0 );
	runtime_assert( factory_ != 0 );

	using core::kinematics::MoveMapOP;
	using core::kinematics::MoveMap;

	// non-initialization failsafe
	if ( ! minmover_->movemap() ) minmover_->movemap( core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new MoveMap ) ) );

	//clone the MinMover's MoveMap
	MoveMapOP mm( new MoveMap( *( minmover_->movemap() ) ) );
	MoveMapOP mm_copy( new MoveMap( *( minmover_->movemap() ) ) );

	//generate task
	using core::pack::task::PackerTaskOP;
	PackerTaskOP task( factory_->create_task_and_apply_taskoperations( pose ) );

	//modify movemap by task
	//  core::kinematics::modify_movemap_from_packertask( *mm, *task );
	Size const nres( task->total_residue() );

	mm->set_jump( jump_ );
	for ( Size i(1); i <= nres; ++i ) {
		if ( task->design_residue( i ) ) {
			// the MoveMap initializes to false for all degrees of freedom
			// this class only turns on minimization at packable dofs, it does not turn them off
			if ( chi_ ) mm->set_chi( i, chi_ );
			if ( bb_ ) mm->set_bb( i, bb_ );
		} else if ( task->pack_residue( i ) ) {
			if ( chi_ ) mm->set_chi( i, chi_ );
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

protocols::moves::MoverOP TaskAwareMinMover::fresh_instance() const { return protocols::moves::MoverOP( new TaskAwareMinMover ); }
protocols::moves::MoverOP TaskAwareMinMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::TaskAwareMinMover( *this ) ); }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
TaskAwareMinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
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
	if ( tag->hasOption( "jump" ) ) {
		std::string const jump( tag->getOption< std::string >( "jump" ) );
		if ( jump != "0" && jump != "1" ) {
			utility_exit_with_message( "TaskAwareMinMover only knows how to interpret jump=1(all jumps true) or jump=0 (false). I got jump = "+jump );
		}
	}
	jump_ = tag->getOption< bool >( "jump", false );
	minmover_ = protocols::simple_moves::MinMoverOP( new MinMover );
	// call to MinMover::parse_my_tag avoided here to prevent collision of chi and bb tag options
	minmover_->parse_opts( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
	minmover_->score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap) );
}

void TaskAwareMinMover::parse_def( utility::lua::LuaObject const & def,
	utility::lua::LuaObject const & score_fxns,
	utility::lua::LuaObject const & tasks,
	protocols::moves::MoverCacheSP cache ) {
	if ( def["chi"] ) {
		chi_ = def["chi"].to<bool>();
	}
	if ( def["bb"] ) {
		bb_ = def["bb"].to<bool>();
	}
	jump_ = def["jump"] ? def["jump"].to<bool>() : false;
	minmover_ = protocols::simple_moves::MinMoverOP( new MinMover );
	// interesting how this doesn't read a movemap....
	minmover_->parse_def_opts( def, score_fxns, tasks, cache );
	if ( def["tasks"] ) {
		core::pack::task::TaskFactoryOP new_task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
		if ( new_task_factory == 0 ) return;
		factory_ =  new_task_factory;
	}
	if ( def["scorefxn"] ) {
		minmover_->score_function( protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns ) );
	} else {
		minmover_->score_function( score_fxns["score12"].to<core::scoring::ScoreFunctionSP>()->clone()  );
	}
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
TaskAwareMinMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	using namespace core::pack::task;
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0 ) return;
	factory_ = new_task_factory;
}

}//moves
}//protocols

