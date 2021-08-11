// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/VirtualResidueSelector.fwd.hh
/// @brief a residue selector that wraps is_virtual_residue
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_VirtualResidueSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_VirtualResidueSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace select {
namespace residue_selector {

class VirtualResidueSelector;

typedef utility::pointer::shared_ptr< VirtualResidueSelector > VirtualResidueSelectorOP;
typedef utility::pointer::shared_ptr< VirtualResidueSelector const > VirtualResidueSelectorCOP;

} //protocols
} //pose_sewing
} //residue_selectors

#endif //INCLUDED_core_select_residue_selector_VirtualResidueSelector_fwd_hh
