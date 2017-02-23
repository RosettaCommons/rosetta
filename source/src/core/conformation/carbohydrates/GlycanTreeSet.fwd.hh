// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core.pose.carbohydrates.GlycanTreeSet.fwd.hh
/// @brief Class to store info on all glycan trees of a pose
/// @author raemisch (raemisch@scripps.edu)


#ifndef INCLUDED_core_conformation_carbohydrates_GlycanTreeSet_fwd_hh
#define INCLUDED_core_conformation_carbohydrates_GlycanTreeSet_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace core {
namespace conformation {
namespace carbohydrates {

class GlycanTreeSet;

typedef utility::pointer::shared_ptr< GlycanTreeSet > GlycanTreeSetOP;
typedef utility::pointer::shared_ptr< GlycanTreeSet const > GlycanTreeSetCOP;



} //core
} //chemical
} //carbohydrates


#endif //INCLUDED_core_pose_carbohydrates_GlycanTreeSet_fwd_hh





