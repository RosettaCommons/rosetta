// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

static basic::Tracer TR( "core.pack.task.operation.TaskOperation" );

TaskOperation::~TaskOperation() = default;

void TaskOperation::parse_tag( TagCOP tag,  DataMap & )
{
	TR << "TaskOperation::parse_tag method called with no effect";
	if ( tag.get() != nullptr ) TR << " for Tag with type " << tag->getName();
	TR << ". Probably due to (un/mis)implemented virtual method in derived class." << std::endl;
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
