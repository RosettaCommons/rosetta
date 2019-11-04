// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/RestrictInteractionGraphThreadsOperation.cc
/// @brief A task operation that restricts the number of threads allowed for interaction graph computation.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <core/pack/task/operation/RestrictInteractionGraphThreadsOperation.hh>
#include <core/pack/task/operation/RestrictInteractionGraphThreadsOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

static basic::Tracer TR( "core.pack.task.operation.RestrictInteractionGraphThreadsOperation" );

namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief Default constructor.
RestrictInteractionGraphThreadsOperation::RestrictInteractionGraphThreadsOperation() = default;

/// @brief Initialization constructor.
/// @param[in] thread_limit_in The number of threads to limit the packer to using. If this is set to zero,
/// this operation does nothing.  Otherwise, it *reduces* the allowed number of threads to the specified value.
/// If the allowed number of threads is less than the specified value, it does nothing.
RestrictInteractionGraphThreadsOperation::RestrictInteractionGraphThreadsOperation(
	core::Size const thread_limit_in
) :
	TaskOperation(),
	num_threads_(thread_limit_in)
{}

/// @brief Default copy constructor.
RestrictInteractionGraphThreadsOperation::RestrictInteractionGraphThreadsOperation( RestrictInteractionGraphThreadsOperation const & /*src*/ ) = default;

/// @brief Destructor.
RestrictInteractionGraphThreadsOperation::~RestrictInteractionGraphThreadsOperation() {}

/// @brief Clone operation: make a copy of this object, and return a smart pointer to the copy.
TaskOperationOP
RestrictInteractionGraphThreadsOperation::clone() const {
	return TaskOperationOP( utility::pointer::make_shared< RestrictInteractionGraphThreadsOperation >( *this ) );
}

/// @brief Configure from a RosettaScripts XML tag.
void
RestrictInteractionGraphThreadsOperation::parse_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
){
	set_thread_limit( tag->getOption<core::Size>( "thread_limit", thread_limit() ) );
}

/// @brief Alter a PackerTask by reducing the number of allowed threads for packer setup to the number
/// specified by this TaskOperation.
/// @details Does nothing if this TaskOperation's allowed threads are set to zero or if the number already
/// allowed in the PackerTask is less than the number allowed by the TaskOperation.
void
RestrictInteractionGraphThreadsOperation::apply(
	core::pose::Pose const & ,
	core::pack::task::PackerTask &task
) const {
	if ( thread_limit() != 0 ) {
		task.limit_ig_setup_threads( thread_limit() );
	}
}

/// @brief Return the name used to construct this TaskOperation from an XML file.
std::string
RestrictInteractionGraphThreadsOperation::keyname() {
	return "RestrictInteractionGraphThreadsOperation";
}

/// @brief Describe the format of XML file used to initialize this TaskOperation.
void
RestrictInteractionGraphThreadsOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"thread_limit",
		xsct_non_negative_integer,
		"The maximum number of threads allowed for interaction graph precomputation.  Setting this to 0 imposes no limit.  Note that this task operation will only lower the number of allowed threads.  That is, if other flags or task operations limit the allowed threads to fewer than this task operation allows, the minimum value will be used.",
		"0"
	);

	task_op_schema_w_attributes( xsd, keyname(), attlist, "The RestrictInteractionGraphThreadsOperation limits the number of threads allowed for interaction graph pre-computation.  By default, all available threads are used.  Note that this only has an effect in the multi-threaded build of Rosetta (built with the \"extras=cxx11thread\" option).\nThis TaskOperation was written on Saturday, 19 October 2019 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." );
}

/// @brief Set the maximum number of threads that this TaskOperation will allow for packer setup.
/// @details Set this to zero to indicate no limit.
void
RestrictInteractionGraphThreadsOperation::set_thread_limit(
	core::Size const setting
) {
#ifndef MULTI_THREADED
	runtime_assert_string_msg( setting < 2, "Error in RestrictInteractionGraphThreadsOperation::set_thread_limit(): In the non-threaded build of Rosetta, the thread limit must be 0 or 1.  Please compile with the \"extras=cxx11thread\" option to enable multi-threading." );
#endif
	TR << "Configuring RestrictInteractionGraphThreadsOperation to allow " << setting << " threads for interaction graph computation." << std::endl;
	num_threads_ = setting;
}

TaskOperationOP
RestrictInteractionGraphThreadsOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< RestrictInteractionGraphThreadsOperation >();
}

std::string
RestrictInteractionGraphThreadsOperationCreator::keyname() const
{
	return RestrictInteractionGraphThreadsOperation::keyname();
}

void
RestrictInteractionGraphThreadsOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictInteractionGraphThreadsOperation::provide_xml_schema( xsd );
}


} //operation
} //task
} //pack
} //core
