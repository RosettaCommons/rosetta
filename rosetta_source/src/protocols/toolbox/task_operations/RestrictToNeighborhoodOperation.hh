// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh
/// @brief  TaskOperation class that finds a neighborhood and makes it mobile in the PackerTask
/// @author Steven Lewis smlewi@unc.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToNeighborhoodOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToNeighborhoodOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near a neighborhood.
class RestrictToNeighborhoodOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;
	typedef std::set< core::Size > SizeSet;

	RestrictToNeighborhoodOperation();

	RestrictToNeighborhoodOperation( std::set< core::Size > const & central_residues );

	RestrictToNeighborhoodOperation( std::string const & calculator );

	virtual ~RestrictToNeighborhoodOperation();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

private:
	///@brief constructor helper function - makes the PoseMetricCalculator
	void make_calculator( std::set< core::Size > const & central_residues );

	///@brief constructor helper function - names the PoseMetricCalculator
	void make_name( std::set< core::Size > const & central_residues );

	std::string calculator_name_;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictToNeighborhoodOperation_HH
