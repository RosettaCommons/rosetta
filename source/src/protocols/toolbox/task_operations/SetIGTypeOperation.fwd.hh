// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SetIGTypeOperation.fwd.hh
///
/// @brief Task operation to set Interaction graph Type
/// @author Sagar Khare



#ifndef INCLUDED_protocols_toolbox_task_operations_SetIGTypeOperation_FWD_HH
#define INCLUDED_protocols_toolbox_task_operations_SetIGTypeOperation_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class SetIGTypeOperation;
typedef utility::pointer::shared_ptr< SetIGTypeOperation > SetIGTypeOperationOP;
typedef utility::pointer::shared_ptr< SetIGTypeOperation const > SetIGTypeOperationCOP;

} //task_operations
} //toolbox
} //protocols

#endif


