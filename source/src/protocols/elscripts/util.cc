// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/elscripts/util.hh
/// @brief Utility functions useful in elscripts
/// @author ken Jung

// Unit Headers
#include <protocols/elscripts/util.hh>

#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/exit.hh>

static thread_local basic::Tracer TR( "protocols.elscripts.util" );

namespace protocols {
namespace elscripts {

using namespace utility;
using namespace utility::lua;

core::scoring::ScoreFunctionOP
parse_scoredef( LuaObject const & scoredef,
		LuaObject const & score_fxns ) {

	std::string scorefxn_name = scoredef.to<std::string>();
	for (LuaIterator i=score_fxns.begin(), end; i != end; ++i) {
		if( i.skey() == scorefxn_name ) {
			return (*i).to<core::scoring::ScoreFunctionSP>()->clone();
			break;
		}
	}
	utility_exit_with_message("ScoreFunction " + scorefxn_name + " not found in score_fxns map.");
}

core::pack::task::TaskFactoryOP
parse_taskdef( LuaObject const & taskdef,
		LuaObject const & tasks ) {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

  TaskFactoryOP new_task_factory( new TaskFactory );
	for (LuaIterator i=taskdef.begin(), end; i != end; ++i) {
		std::string task_name = (*i).to<std::string>();
		bool taskfound = false;
		for (LuaIterator j=tasks.begin(), end; j != end; ++j) {
			if( j.skey() == task_name ) {
				new_task_factory->push_back( (*j).to<TaskOperationOP>() );
				taskfound = true;
				break;
			}
		}
		if( !taskfound ) utility_exit_with_message("TaskOperation " + task_name + " not found in tasks map.");
	}
  return new_task_factory;
}

core::pack::task::TaskFactorySP
sp_parse_taskdef( LuaObject const & taskdef,
		LuaObject const & tasks ) {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	TaskFactoryOP tmpop = parse_taskdef( taskdef, tasks );
	TaskFactorySP tmpsp( tmpop.get() );
	tmpop.reset(); // No relinquish_ownership in std::shared_ptr
	return tmpsp;
}

void parse_movemapdef( LuaObject const & movemapdef, core::kinematics::MoveMapOP mm ) {
	for (LuaIterator i=movemapdef.begin(), end; i != end; ++i) {
		switch( (*i).size() ) {
			case 2:
				{
					core::Size const num( (*i)[1].to<core::Size>() );
					bool const setting( (*i)[2].to<bool>() );
					num == 0 ? mm->set_jump( setting ) : mm->set_jump( num, setting );
					break;
				}
			case 4:
				{
					core::Size const begin( (*i)[1].to<core::Size>() );
					core::Size const end( (*i)[2].to<core::Size>() );
					runtime_assert( end >= begin );
					bool const chi( (*i)["chi"].to<bool>() );
					bool const bb( (*i)["bb"].to<bool>() );
					for( core::Size i( begin ); i <= end; ++i ){
						mm->set_chi( i, chi );
						mm->set_bb( i, bb );
					}
					break;
				}
		}
	}
}

} //elscripts
} //protocols
