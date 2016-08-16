// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNode.fwd.hh
/// @brief  Forward declaration for inverse rotamer target
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNode_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNode_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class InvrotTreeNode;

typedef utility::pointer::shared_ptr< InvrotTreeNode > InvrotTreeNodeOP;
typedef utility::pointer::shared_ptr< InvrotTreeNode const > InvrotTreeNodeCOP;

}
}
}

#endif
