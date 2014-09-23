// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictByCalculatorsOperation.cc
/// @brief  A class that applies arbitrary calculators (whose calculations return std::set< core::Size >) to restrict a PackerTask
/// @author Steven Lewis smlewi@unc.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>

// AUTO-REMOVED #include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


// C++ Headers

using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictByCalculatorsOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

RestrictByCalculatorsOperation::RestrictByCalculatorsOperation() {}

RestrictByCalculatorsOperation::RestrictByCalculatorsOperation(utility::vector1< calc_calcn > const & calcs_and_calcns)
	: parent(), calcs_and_calcns_(calcs_and_calcns)
{
	//I suppose you could reasonably create this object BEFORE the calculator was generated/registered
// 	for( core::Size i(1); i <= calcs_and_calcns_.size(); ++i){
// 		if( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calcs_and_calcns_[i].first ) ){
// 			utility_exit_with_message("In RestrictByCalculatorsOperation, calculator " + calcs_and_calcns_[i].first + " does not exist.");
// 		}
// 	}
}

RestrictByCalculatorsOperation::~RestrictByCalculatorsOperation() {}

core::pack::task::operation::TaskOperationOP
RestrictByCalculatorsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictByCalculatorsOperation );
}

///@details be warned if you use clone that you'll not get new calculators
core::pack::task::operation::TaskOperationOP RestrictByCalculatorsOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictByCalculatorsOperation( *this ) );
}

void
RestrictByCalculatorsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
 	//vector for filling packertask
 	utility::vector1_bool repack(pose.total_residue(), false);

 	for( core::Size i(1); i <= calcs_and_calcns_.size(); ++i)
		run_calculator(pose, calcs_and_calcns_[i].first, calcs_and_calcns_[i].second, repack);

 	task.restrict_to_residues(repack);
	return;
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
