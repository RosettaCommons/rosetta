// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperation.cc
/// @brief  An operation to perform on a packer task --
///         usually, by a PackerTaskFactory right after the task's construction
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// only generalized base classes go here. TaskOperations that actually do things do not belong here.

// Unit Headers
#include <core/pack/task/operation/TaskOperation.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

#ifdef USELUA
void lregister_TaskOperation( lua_State * lstate ) {
	luabind::module(lstate, "core")
	[
		luabind::namespace_("pack")
		[
			luabind::namespace_("task")
			[
				luabind::namespace_("operation")
				[
					luabind::class_<TaskOperation>("TaskOperation")
				]
			]
		]
	];
}
#endif

static THREAD_LOCAL basic::Tracer TR( "core.pack.task.operation.TaskOperation" );

TaskOperation::~TaskOperation() {}

void TaskOperation::parse_tag( TagCOP tag,  DataMap & )
{
	TR << "TaskOperation::parse_tag method called with no effect";
	if ( tag.get() != NULL ) TR << " for Tag with type " << tag->getName();
	TR << ". Probably due to (un/mis)implemented virtual method in derived class." << std::endl;
}

void TaskOperation::parse_def( utility::lua::LuaObject const & /*def*/ ) {
	utility_exit_with_message("This TaskOperation has not implemented parse_def()");
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
