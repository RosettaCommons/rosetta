// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Declaration of the base class for ResLvlTaskOperation factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperationRegistrator_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperationRegistrator_hh

// Package headers
#include <utility/factory/WidgetRegistrator.hh>

#include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief This templated class will register an instance of an
/// ResLvlTaskOperationCreator (class T) with the ResLvlTaskOperationFactory.  It will ensure
/// that no ResLvlTaskOperation creator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ResLvlTaskOperationRegistrator : public utility::factory::WidgetRegistrator< ResLvlTaskOperationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResLvlTaskOperationFactory, T > parent;
public:
	ResLvlTaskOperationRegistrator() : parent() {}
};

}
}
}
}

#endif
