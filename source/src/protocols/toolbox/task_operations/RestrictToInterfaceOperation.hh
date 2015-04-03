// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh
/// @brief  TaskOperation class that finds an interface and makes it mobile in the PackerTask
/// @author Steven Lewis smlewi@unc.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToInterfaceOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToInterfaceOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictToInterfaceOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	RestrictToInterfaceOperation( core::Size upper_chain = 1, core::Size lower_chain = 2 );

	RestrictToInterfaceOperation( std::string const & calculator );

	virtual ~RestrictToInterfaceOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

private:
	/// @brief constructor helper function - makes the PoseMetricCalculator
	void make_calculator( core::Size upper_chain, core::Size lower_chain );

	/// @brief constructor helper function - names the PoseMetricCalculator
	void make_name( core::Size upper_chain, core::Size lower_chain );

	std::string calculator_name_;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictToInterfaceOperation_HH
