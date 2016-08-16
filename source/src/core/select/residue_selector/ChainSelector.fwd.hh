// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ChainSelector.fwd.hh
/// @brief  The ChainSelector selects all the residues from a given chain, given either as a number
///         or a character.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ChainSelector_FWD_HH
#define INCLUDED_core_select_residue_selector_ChainSelector_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace residue_selector {

class ChainSelector;

typedef utility::pointer::shared_ptr< ChainSelector > ChainSelectorOP;
typedef utility::pointer::shared_ptr< ChainSelector const > ChainSelectorCOP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
