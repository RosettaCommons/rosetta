// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/DesignAroundOperation.fwd.hh
/// @brief
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_task_operations_DesignAroundOperation_fwd_hh
#define INCLUDED_protocols_task_operations_DesignAroundOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace task_operations {

class DesignAroundOperation;

typedef utility::pointer::shared_ptr< DesignAroundOperation > DesignAroundOperationOP;
typedef utility::pointer::shared_ptr< DesignAroundOperation const > DesignAroundOperationCOP;

} //namespace protocols
} //namespace task_operations

#endif // INCLUDED_protocols_TaskOperations_DesignAroundOperation_FWD_HH

