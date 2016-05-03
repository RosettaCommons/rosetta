// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/pack/task/operation/ReplicateTask.cc
/// @details Replicates logic from one pose onto another.  Poses must be the same size!
/// Only copies logic for being_packed and being_designed.
/// Use this to retain task information from round to round and keep
/// sequence information from being corrupted in the task.
/// @author stranges

#include <core/pack/task/operation/ReplicateTask.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>

#include <utility/exit.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

/// @details empty constructor need to call set_native_task to make it work
ReplicateTask::ReplicateTask() : parent(){}

/// @details constructor using an established PackerTask
ReplicateTask::ReplicateTask( core::pack::task::PackerTaskOP native_task ) :
	parent(), native_task_( native_task->clone() ){}

/// @details contstructor that uses an TaskFactory and applies it to the native to get the task
ReplicateTask::ReplicateTask(core::pose::Pose & native_pose, core::pack::task::TaskFactoryOP task_factory ) :
	parent()
{
	//convert task factory into task for general use
	native_task_ =  task_factory->create_task_and_apply_taskoperations( native_pose ) ;
}

ReplicateTask::~ReplicateTask(){}

task::operation::TaskOperationOP ReplicateTask::clone() const
{
	return task::operation::TaskOperationOP( new ReplicateTask( *this ) );
}

void
ReplicateTask::apply(
	pose::Pose const & pose,
	task::PackerTask & task) const
{
	// note:  poses must be the same size as this copies the task logic
	// on a per residue basis
	runtime_assert( pose.total_residue() == native_task_->total_residue() );
	//for all the residues in the pose copy the basic task logic
	for ( Size ii = 1; ii <= pose.total_residue(); ii++ ) {
		//if not being designed then restrict to repacking
		if ( !native_task_->nonconst_residue_task( ii ).being_designed() ) {
			task.nonconst_residue_task( ii ).restrict_to_repacking();
		}
		//task.nonconst_residue_task( ii ).add_behavior( "NATAA" );

		//if not being packed at all then prevent from repacking
		if ( !native_task_->nonconst_residue_task( ii ).being_packed() ) {
			task.nonconst_residue_task( ii ).prevent_repacking();
		}
		//task.nonconst_residue_task( ii ).add_behavior( "NATRO" );

	} //end loop over all residues
} //end apply

void
ReplicateTask::set_native_task( core::pack::task::PackerTaskOP native_task){
	native_task_ = native_task;
}

//private:

} //operation
} //pack
} //task
} //core
