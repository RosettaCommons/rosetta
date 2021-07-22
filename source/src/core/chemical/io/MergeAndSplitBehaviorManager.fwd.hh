// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/io/MergeAndSplitBehaviorManager.fwd.hh
/// @brief   Forward declarations for MergeAndSplitBehaviorManager.


#ifndef INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_FWD_HH
#define INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace io {

/// @brief  A class for storing instructions on how to merge or split residues in a particular TypeSet when loading
///         them from the PDB.
class MergeAndSplitBehaviorManager;

typedef utility::pointer::shared_ptr< MergeAndSplitBehaviorManager > MergeAndSplitBehaviorManagerOP;
typedef utility::pointer::shared_ptr< MergeAndSplitBehaviorManager const > MergeAndSplitBehaviorManagerCOP;

}  // namespace io
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_FWD_HH
