// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/DesignRestrictions.fwd.hh
/// @brief Combine a set of ResidueSelectors and ResidueLevelTaskOperations concisely, using boolean logic to join selectors concisely
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_task_operation_DesignRestrictions_fwd_hh
#define INCLUDED_core_pack_task_operation_DesignRestrictions_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace task {
namespace operation {

class DesignRestrictions;

typedef utility::pointer::shared_ptr< DesignRestrictions > DesignRestrictionsOP;
typedef utility::pointer::shared_ptr< DesignRestrictions const > DesignRestrictionsCOP;

} //core
} //pack
} //task
} //operation

#endif //INCLUDED_core_pack_task_operation_DesignRestrictions_fwd_hh
