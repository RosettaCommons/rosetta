// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.fwd.hh
/// @brief Restrict loop residue identities to those commonly found in nature
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperation_fwd_hh
#define INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperation_fwd_hh

// utilty headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace denovo_design {
namespace task_operations {

// Forward declaration
class ConsensusLoopDesignOperation;

// Types
typedef utility::pointer::shared_ptr< ConsensusLoopDesignOperation > ConsensusLoopDesignOperationOP;
typedef utility::pointer::shared_ptr< ConsensusLoopDesignOperation const > ConsensusLoopDesignOperationCOP;

}
}
}
#endif


