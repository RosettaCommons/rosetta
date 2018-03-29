// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictResiduesToRepackingOperation.fwd.hh
/// @brief
/// @author Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_protocols_task_operations_RestrictResiduesToRepackingOperation_fwd_hh
#define INCLUDED_protocols_task_operations_RestrictResiduesToRepackingOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace task_operations {

class RestrictResiduesToRepackingOperation;

typedef utility::pointer::shared_ptr< RestrictResiduesToRepackingOperation > RestrictResiduesToRepackingOperationOP;
typedef utility::pointer::shared_ptr< RestrictResiduesToRepackingOperation const > RestrictResiduesToRepackingOperationCOP;

} //namespace protocols
} //namespace task_operations

#endif // INCLUDED_protocols_TaskOperations_RestrictResiduesToRepackingOperation_FWD_HH

