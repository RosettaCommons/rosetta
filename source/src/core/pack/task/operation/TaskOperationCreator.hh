// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/methods/TaskOperationCreator.hh
/// @brief  Declaration of the base class for TaskOperation factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_TaskOperationCreator_hh
#define INCLUDED_core_pack_task_operation_TaskOperationCreator_hh

// Unit headers
#include <core/pack/task/operation/TaskOperationCreator.fwd.hh>

// Package headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief The TaskOperationCreator class's responsibilities are to create
/// on demand a new TaskOperation class.
/// The TaskOperationCreator must register itself with the TaskOperationFactory
/// at load time (before main() begins) so that the TaskOperationFactory is ready
/// to start creating TaskOperations by the time any protocol
/// requests one.
class TaskOperationCreator : public utility::pointer::ReferenceCount
{
public:
	virtual ~TaskOperationCreator() {}

	/// @brief Instantiate a new TaskOperation
	virtual TaskOperationOP create_task_operation() const = 0;
	virtual std::string keyname() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

}
}
}
}

#endif
