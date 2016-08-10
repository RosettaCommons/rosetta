// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.fwd.hh
/// @brief Selects residues that are involved in strand-strand pairings
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_residue_selectors_PairedSheetResidueSelector_fwd_hh
#define INCLUDED_protocols_denovo_design_residue_selectors_PairedSheetResidueSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace denovo_design {
namespace residue_selectors {

class PairedSheetResidueSelector;

typedef utility::pointer::shared_ptr< PairedSheetResidueSelector > PairedSheetResidueSelectorOP;
typedef utility::pointer::shared_ptr< PairedSheetResidueSelector const > PairedSheetResidueSelectorCOP;

} //protocols
} //denovo_design
} //residue_selectors

#endif //INCLUDED_protocols_denovo_design_residue_selectors_PairedSheetResidueSelector_fwd_hh

