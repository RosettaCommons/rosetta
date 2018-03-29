// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RetrieveStoredTaskOperation.fwd.hh
/// @brief
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_protocols_task_operations_RetrieveStoredTaskOperation_fwd_hh
#define INCLUDED_protocols_task_operations_RetrieveStoredTaskOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace task_operations {

class RetrieveStoredTaskOperation;

typedef utility::pointer::shared_ptr< RetrieveStoredTaskOperation > RetrieveStoredTaskOperationOP;
typedef utility::pointer::shared_ptr< RetrieveStoredTaskOperation const > RetrieveStoredTaskOperationCOP;

} //namespace task_operations
} //namespace protocols

#endif // INCLUDED_protocols_task_operations_RetrieveStoredTaskOperation_FWD_HH

