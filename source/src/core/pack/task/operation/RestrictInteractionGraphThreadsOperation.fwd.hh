// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/RestrictInteractionGraphThreadsOperation.fwd.hh
/// @brief A task operation that restricts the number of threads allowed for interaction graph computation.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_core_pack_task_operation_RestrictInteractionGraphThreadsOperation_fwd_hh
#define INCLUDED_core_pack_task_operation_RestrictInteractionGraphThreadsOperation_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictInteractionGraphThreadsOperation;

using RestrictInteractionGraphThreadsOperationOP = utility::pointer::shared_ptr< RestrictInteractionGraphThreadsOperation >;
using RestrictInteractionGraphThreadsOperationCOP = utility::pointer::shared_ptr< RestrictInteractionGraphThreadsOperation const >;

} //operation
} //task
} //pack
} //core

#endif //INCLUDED_core_pack_task_operation_RestrictInteractionGraphThreadsOperation_fwd_hh
