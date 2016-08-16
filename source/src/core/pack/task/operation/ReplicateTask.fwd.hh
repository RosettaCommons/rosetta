// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ReplicateTask.fwd.hh
/// @brief  forward declaration for ReplicateTask
/// @author Ben Stranges (stranges@unc.edu)

#ifndef INCLUDED_core_pack_task_operation_ReplicateTask_fwd_hh
#define INCLUDED_core_pack_task_operation_ReplicateTask_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief forward declaration for ReplicateTask
class ReplicateTask;


/// @brief ReplicateTask owning pointer
typedef utility::pointer::shared_ptr< ReplicateTask > ReplicateTaskOP;


/// @brief ReplicateTask const owning pointer
typedef utility::pointer::shared_ptr< ReplicateTask const > ReplicateTaskCOP;


/// @brief ReplicateTask owning pointer
typedef utility::pointer::weak_ptr< ReplicateTask > ReplicateTaskAP;


/// @brief ReplicateTask const owning pointer
typedef utility::pointer::weak_ptr< ReplicateTask const > ReplicateTaskCAP;


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_ReplicateTask_FWD_HH */
