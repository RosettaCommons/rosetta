// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NumNeighborsSelector.fwd.hh
/// @brief  Forward declaration of a class that identifies residues with at least X neighbors within a cutoff distance D
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_NumNeighborsSelector_FWD_HH
#define INCLUDED_core_select_residue_selector_NumNeighborsSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace residue_selector {

class NumNeighborsSelector;

typedef utility::pointer::shared_ptr< NumNeighborsSelector > NumNeighborsSelectorOP;
typedef utility::pointer::shared_ptr< NumNeighborsSelector const > NumNeighborsSelectorCOP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
