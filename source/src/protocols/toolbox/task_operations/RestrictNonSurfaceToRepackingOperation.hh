// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh
/// @brief  TaskOperation class that restrics buried and boundary positions to repacking leaving surface positions designable
/// @author Ron Jacak ronj@unc.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictNonSurfaceToRepackingOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictNonSurfaceToRepackingOperation_hh

// Unit Headers
#include <core/pack/task/operation/TaskOperation.hh> // abstract base class
#include <protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to allow design of only surface exposed positions
class RestrictNonSurfaceToRepackingOperation : public core::pack::task::operation::TaskOperation {

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;

public:
	RestrictNonSurfaceToRepackingOperation();
	RestrictNonSurfaceToRepackingOperation( core::Size nb_cutoff );
	virtual ~RestrictNonSurfaceToRepackingOperation();

	virtual TaskOperationOP clone() const;
	void surface_exposed_nb_cutoff( core::Size const nb_count );
	virtual	void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	virtual void parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & dm );

private:
	core::Size surface_exposed_nb_cutoff_; // 16, by default

};


} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictNonSurfaceToRepackingOperation_HH
