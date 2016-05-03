// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/methods/ResLvlTaskOperationCreator.hh
/// @brief  Declaration of the base class for ResLvlTaskOperation factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth


#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperationCreator_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperationCreator_hh

// Unit headers
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>

// Package headers
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief The ResLvlTaskOperationCreator class's responsibilities are to create
/// on demand a new ResLvlTaskOperation class.
/// The ResLvlTaskOperationCreator must register itself with the ResLvlTaskOperationFactory
/// at load time (before main() begins) so that the ResLvlTaskOperationFactory is ready
/// to start creating ResLvlTaskOperations by the time any protocol
/// requests one.
class ResLvlTaskOperationCreator : public utility::pointer::ReferenceCount
{
public:
	virtual ~ResLvlTaskOperationCreator() {}

	/// @brief Instantiate a new ResLvlTaskOperation
	virtual
	ResLvlTaskOperationOP
	create_res_level_task_operation() const = 0;

	virtual std::string keyname() const = 0;

	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

}
}
}
}

#endif
