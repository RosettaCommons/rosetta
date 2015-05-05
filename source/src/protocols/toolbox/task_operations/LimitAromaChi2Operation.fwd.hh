// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/LimitAromaChi2Operation.hh
/// @brief  rotamer set operation forward declaration
/// @author Nobuyasu Koga (nobuyasu@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_LimitAromaChi2Operation_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_LimitAromaChi2Operation_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {


class LimitAromaChi2_RotamerSetOperation;
class LimitAromaChi2Operation;

typedef utility::pointer::shared_ptr< LimitAromaChi2_RotamerSetOperation > LimitAromaChi2_RotamerSetOperationOP;
typedef utility::pointer::shared_ptr< LimitAromaChi2_RotamerSetOperation const > LimitAromaChi2_RotamerSetOperationCOP;

typedef utility::pointer::shared_ptr< LimitAromaChi2Operation > LimitAromaChi2OperationOP;
typedef utility::pointer::shared_ptr< LimitAromaChi2Operation const > LimitAromaChi2OperationCOP;


}
}
}


#endif
