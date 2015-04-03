// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/pack/task/operation/ReplicateTask.hh
/// @details Replicates logic from one pose onto another.  Poses must be the same size!
/// Only copies logic for being_packed and being_designed
/// @author stranges

#ifndef INCLUDED_core_pack_task_operation_ReplicateTask_hh
#define INCLUDED_core_pack_task_operation_ReplicateTask_hh

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/ReplicateTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>


// for parsing
#include <utility/tag/Tag.fwd.hh>
    //#include <utility/vector1.hh>

#include <core/types.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

class ReplicateTask : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
public:
	//constructors
	/// @brief empty constructor need to call set_native_task(task) to make it work
	ReplicateTask();
	/// @brief actual useful constructors
	ReplicateTask( core::pack::task::PackerTaskOP  native_task );
	ReplicateTask( core::pose::Pose & native_pose, core::pack::task::TaskFactoryOP  task_factory );

	virtual ~ReplicateTask();
	virtual TaskOperationOP clone() const;

	virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	/// Does NOT Work! DO NOT USE parse_tag here
	virtual void parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & );
	//void symmetric_task( core::pose::Pose const & pose, task::PackerTask & task ) const;
	void set_native_task( core::pack::task::PackerTaskOP native_task);


private:
	core::pack::task::PackerTaskOP native_task_;

};//end of class ReplicateTask

} //operation
} //pack
} //task
} //core

#endif
