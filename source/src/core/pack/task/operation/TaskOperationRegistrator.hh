// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Declaration of the base class for TaskOperation factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_TaskOperationRegistrator_hh
#define INCLUDED_core_pack_task_operation_TaskOperationRegistrator_hh

// Package headers
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief This templated class will register an instance of an
/// TaskOperationCreator (class T) with the TaskOperationFactory.  It will ensure
/// that no TaskOperation creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class TaskOperationRegistrator : public utility::factory::WidgetRegistrator< TaskOperationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< TaskOperationFactory, T > parent;
public:
	TaskOperationRegistrator() : parent() {}
};

}
}
}
}

#endif
