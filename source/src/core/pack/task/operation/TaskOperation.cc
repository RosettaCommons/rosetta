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

/// @brief Does this task operation provide information about how to cite it?
/// @details Defaults to false.  Derived classes may override this to provide citation info.  If set to
/// true, the provide_citation_info() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
TaskOperation::task_operation_provides_citation_info() const {
	return false;
}

/// @brief Provide the citation.
/// @returns A vector of citation collections.  This allows the task operation to provide citations for
/// itself and for any modules that it invokes.
/// @details The default implementation of this function provides an empty vector.  It may be
/// overriden by task operations wishing to provide citation information.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
TaskOperation::provide_citation_info() const {
	return utility::vector1< basic::citation_manager::CitationCollectionCOP >();
}

/// @brief Does this task_operation indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @details Defaults to false.  Derived classes may override this to provide authorship info.  If set to
/// true, the provide_authorship_info_for_unpublished() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
TaskOperation::task_operation_is_unpublished() const {
	return false;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @returns A list of pairs of (author, e-mail address).  Empty list if not unpublished.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
TaskOperation::provide_authorship_info_for_unpublished() const {
	return utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >();
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
