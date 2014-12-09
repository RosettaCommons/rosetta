// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SetIGTypeOperation.cc
///
/// @brief Task operation to set interaction graph type (linmem, lazy or double lazy)
/// @author Sagar Khare

//Unit Headers
#include <protocols/toolbox/task_operations/SetIGTypeOperation.hh>
#include <protocols/toolbox/task_operations/SetIGTypeOperationCreator.hh>

//Core Headers
#include <core/pack/task/PackerTask.hh>

//Utility Headers
#include <utility/tag/Tag.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

SetIGTypeOperation::SetIGTypeOperation():
	lin_mem_(false),
	lazy_(false),
	double_lazy_(false)
{}

SetIGTypeOperation::~SetIGTypeOperation(){}

core::pack::task::operation::TaskOperationOP SetIGTypeOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new SetIGTypeOperation( *this ) );
}

void SetIGTypeOperation::apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const
{
	if (lin_mem_) task.or_linmem_ig( true );
	else if (lazy_) task.or_lazy_ig( true );
	else if (double_lazy_) task.or_double_lazy_ig( true );
}

void SetIGTypeOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	lin_mem_ = tag->getOption< bool >("lin_mem_ig", false);
	lazy_ = tag->getOption< bool >("lazy_ig", false);
	double_lazy_ = tag->getOption< bool >("double_lazy_ig", false);
}

core::pack::task::operation::TaskOperationOP SetIGTypeOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SetIGTypeOperation );
}

} //task_operations
} //toolbox
} //protocols
