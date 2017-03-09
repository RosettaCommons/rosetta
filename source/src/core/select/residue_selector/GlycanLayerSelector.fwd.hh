// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/GlycanLayerSelector.fwd.hh
/// @brief A selector for choosing glycan residues based on their layer - as measured by the residue distance to the start of the glycan tree.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_GlycanLayerSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_GlycanLayerSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace select {
namespace residue_selector {

class GlycanLayerSelector;

typedef utility::pointer::shared_ptr< GlycanLayerSelector > GlycanLayerSelectorOP;
typedef utility::pointer::shared_ptr< GlycanLayerSelector const > GlycanLayerSelectorCOP;

} //core
} //select
} //residue_selector

#endif //INCLUDED_core_select_residue_selector_GlycanLayerSelector_fwd_hh
