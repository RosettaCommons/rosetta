// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanNode.fwd.hh
/// @brief Class to store info a node (residue) within a glycan tree
/// @author raemisch (raemisch@scripps.edu)


#ifndef INCLUDED_core_conformation_carbohydrates_GlycanNode_fwd_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanNode_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace core {
namespace conformation {
namespace carbohydrates {

class GlycanNode;

typedef utility::pointer::shared_ptr< GlycanNode > GlycanNodeOP;
typedef utility::pointer::shared_ptr< GlycanNode const > GlycanNodeCOP;



} //core
} //chemical
} //carbohydrates


#endif //INCLUDED_core_pose_carbohydrates_GlycanNode_fwd_hh





