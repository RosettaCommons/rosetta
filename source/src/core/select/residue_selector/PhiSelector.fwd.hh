// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/PhiSelector.fwd.hh
/// @brief A ResidueSelector that selects alpha-amino acids that are either in the positive
/// phi or negative phi region of Ramachandran space (depending on user preferences).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_core_select_residue_selector_PhiSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_PhiSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace core {
namespace select {
namespace residue_selector {

class PhiSelector;

typedef utility::pointer::shared_ptr< PhiSelector > PhiSelectorOP;
typedef utility::pointer::shared_ptr< PhiSelector const > PhiSelectorCOP;

} //core
} //select
} //residue_selector


#endif //INCLUDED_core_select_residue_selector_PhiSelector_fwd_hh





