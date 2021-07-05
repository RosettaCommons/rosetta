// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RetrieveStoredTaskOperation.cc
/// @brief  Retrieves a stored task from the pose's CacheableData. Must be used in
///         conjunction with the StoredTaskMover. Allows storage/retrieval of a task
///         so that particular sets of residues can be stably addressed throughout
///         the entirety of a RosettaScripts protocol.
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <protocols/task_operations/RetrieveStoredTaskOperation.hh>
#include <protocols/task_operations/RetrieveStoredTaskOperationCreator.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/task_operations/STMStoredTask.hh>

// C++ Headers

//Auto Headers
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


static basic::Tracer TR( "protocols.task_operations.RetrieveStoredTaskOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
RetrieveStoredTaskOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< RetrieveStoredTaskOperation >();
}

void RetrieveStoredTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RetrieveStoredTaskOperation::provide_xml_schema( xsd );
}

std::string RetrieveStoredTaskOperationCreator::keyname() const
{
	return RetrieveStoredTaskOperation::keyname();
}

// @brief default constructor
RetrieveStoredTaskOperation::RetrieveStoredTaskOperation() = default;

// @brief destructor
RetrieveStoredTaskOperation::~RetrieveStoredTaskOperation() = default;

// @brief copy constructor
core::pack::task::operation::TaskOperationOP RetrieveStoredTaskOperation::clone() const
{
	return utility::pointer::make_shared< RetrieveStoredTaskOperation >( *this );
}

// @brief apply function
void
RetrieveStoredTaskOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		utility_exit_with_message("Your pose does not have CacheableData of type STM_STORED_TASKS");
	} else {
		protocols::task_operations::STMStoredTask const & stored_tasks = *( utility::pointer::static_pointer_cast< protocols::task_operations::STMStoredTask const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) );
		if ( !stored_tasks.has_task(task_name_) ) {
			utility_exit_with_message("No stored task with the name " + task_name_ + " found");
		} else {
			task.update_commutative( *( stored_tasks.get_task( task_name_ ) ) );
		}
	}
}

// @brief parse xml
void
RetrieveStoredTaskOperation::parse_tag( TagCOP tag , DataMap & )
{
	task_name_ = tag->getOption< std::string >( "task_name" ) ;
}

void RetrieveStoredTaskOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "task_name", xsct_pose_cached_task_operation, "A previously-stored task operation to retrieve from the datacache of a Pose object." );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Retrieves a task operation that has been cached in the datacache of a Pose object." );
}


} //namespace task_operations
} //namespace protocols
