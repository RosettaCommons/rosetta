// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/PruneBuriedUnsatsOperation.fwd.hh
/// @brief  Forward declaration for PruneBuriedUnsatsOperation
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_protocols_task_operations_PruneBuriedUnsatsOperation_fwd_hh
#define INCLUDED_protocols_task_operations_PruneBuriedUnsatsOperation_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace task_operations {


class PruneBuriedUnsats_RotamerSetsOperation;
class PruneBuriedUnsatsOperation;

typedef utility::pointer::shared_ptr< PruneBuriedUnsats_RotamerSetsOperation > PruneBuriedUnsats_RotamerSetsOperationOP;
typedef utility::pointer::shared_ptr< PruneBuriedUnsats_RotamerSetsOperation const > PruneBuriedUnsats_RotamerSetsOperationCOP;

typedef utility::pointer::shared_ptr< PruneBuriedUnsatsOperation > PruneBuriedUnsatsOperationOP;
typedef utility::pointer::shared_ptr< PruneBuriedUnsatsOperation const > PruneBuriedUnsatsOperationCOP;


}
}


#endif
