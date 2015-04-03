// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictOperationsBase.hh
/// @brief  Base class for PoseMetricCalculator-using TaskOperations
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictOperationsBase_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictOperationsBase_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictOperationsBase.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details This base class defines an interface for TaskOperations which use PoseMetricCalculators to pick out certain residues in a pose.  The apply function takes a vector of PoseMetricCalculators and a vector of calculation queries, and uses those queries to shape the PackerTask.  This was designed to work with NeighborsByDistanceCalculator, NeighborhoodByDistanceCalculator, and InterfaceNeighborDefinitionCalculator; in general it works when the calculators can return MetricValue< std::set< core::Size > > (a std::set of resids).
class RestrictOperationsBase : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;

	RestrictOperationsBase();

	virtual ~RestrictOperationsBase();

	virtual TaskOperationOP clone() const = 0;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const = 0;

protected:
	/// @brief this is the only real function - it takes a calculator name and calculation, and a PackerTask-compatible vector, and flips booleans in the vector according to the calculator
	void
	run_calculator(
		core::pose::Pose const & pose,
		std::string const & calculator,
		std::string const & calculation,
		utility::vector1_bool & residues ) const;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictOperationsBase_HH
