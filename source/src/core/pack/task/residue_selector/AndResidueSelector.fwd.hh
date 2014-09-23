// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/AndResidueSelector.fwd.hh
/// @brief  Forward declaration of a class that combines ResidueSelector logic 
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_AndResidueSelector_FWD_HH
#define INCLUDED_core_pack_task_residue_selector_AndResidueSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

// Project Headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class AndResidueSelector;

typedef utility::pointer::shared_ptr< AndResidueSelector > AndResidueSelectorOP;
typedef utility::pointer::shared_ptr< AndResidueSelector const > AndResidueSelectorCOP;

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
