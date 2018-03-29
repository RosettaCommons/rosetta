// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterfaceTaskOperation.hh
/// @brief Base class for TaskOperations that will work with docking to define an interface
/// @author Brian Weitzner (brian.weitzner@jhu.edu)

#ifndef INCLUDED_protocols_task_operations_InterfaceTaskOperation_hh
#define INCLUDED_protocols_task_operations_InterfaceTaskOperation_hh

#include <protocols/simple_task_operations/InterfaceTaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <core/types.hh>

namespace protocols {
namespace simple_task_operations {


class InterfaceTaskOperation : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
public:
	InterfaceTaskOperation();
	//Copy constructor
	InterfaceTaskOperation( InterfaceTaskOperation const & old_instance );

	virtual ~InterfaceTaskOperation();

	void set_movable_jumps( utility::vector1_int const & movable_jumps );
	void add_movable_jump( int const additional_jump );
	utility::vector1_int const & movable_jumps() const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "InterfaceTaskOperation"; }

private:
	utility::vector1_int movable_jumps_;
};

}
}

#endif
