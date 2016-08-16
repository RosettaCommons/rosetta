// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_SelectByDeltaScoreOperation_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_SelectByDeltaScoreOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class SelectByDeltaScoreOperation;

typedef utility::pointer::shared_ptr< SelectByDeltaScoreOperation > SelectByDeltaScoreOperationOP;

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif // INCLUDED_protocols_toolbox_task_operations_SelectByDeltaScoreOperation_fwd_hh
