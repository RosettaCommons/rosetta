// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/residue_selectors/LigandMetalContactSelector.fwd.hh
/// @brief This residue selector takes a selector or residue number of a ligand and returns any residues in contact with metal atoms in the ligand.
/// @author Allison Watwood (acw538@msstate.edu)

#ifndef INCLUDED_protocols_residue_selectors_LigandMetalContactSelector_fwd_hh
#define INCLUDED_protocols_residue_selectors_LigandMetalContactSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace residue_selectors {

class LigandMetalContactSelector;

typedef utility::pointer::shared_ptr< LigandMetalContactSelector > LigandMetalContactSelectorOP;
typedef utility::pointer::shared_ptr< LigandMetalContactSelector const > LigandMetalContactSelectorCOP;

} //protocols
} //residue_selectors

#endif //INCLUDED_protocols_residue_selectors_LigandMetalContactSelector_fwd_hh
