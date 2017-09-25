// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/FillAUTOTaskOperation.fwd.hh
/// @brief fills the AUTO behavior for all residues in Task. Useful if a protocol expects AUTO-style resfile, but no resfile present.
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_core_pack_task_operation_FillAUTOTaskOperation_fwd_hh
#define INCLUDED_core_pack_task_operation_FillAUTOTaskOperation_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace task {
namespace operation {

class FillAUTOTaskOperation;

typedef utility::pointer::shared_ptr< FillAUTOTaskOperation > FillAUTOTaskOperationOP;
typedef utility::pointer::shared_ptr< FillAUTOTaskOperation const > FillAUTOTaskOperationCOP;

} //core
} //pack
} //task
} //operation

#endif //INCLUDED_core_pack_task_operation_FillAUTOTaskOperation_fwd_hh
