// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/BlockSelector.fwd.hh
/// @brief selectes a specified continuous block of previously selected residues
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_BlockSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_BlockSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace select {
namespace residue_selector {

class BlockSelector;

using BlockSelectorOP = utility::pointer::shared_ptr< BlockSelector >;
using BlockSelectorCOP = utility::pointer::shared_ptr< BlockSelector const >;

} //residue_selectors
} //pose_sewing
} //protocols

#endif //INCLUDED_core_select_residue_selector_BlockSelector_fwd_hh
