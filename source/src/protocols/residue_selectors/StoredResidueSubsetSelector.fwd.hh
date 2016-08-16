// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/StoredResidueSubsetSelector.fwd.hh
/// @brief  Forward declaration of a class that selects residues using a previously stored residue subset
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_residue_selectors_StoredResidueSubsetSelector_FWD_HH
#define INCLUDED_protocols_residue_selectors_StoredResidueSubsetSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace protocols {
namespace residue_selectors {

class StoredResidueSubsetSelector;
typedef utility::pointer::shared_ptr< StoredResidueSubsetSelector > StoredResidueSubsetSelectorOP;
typedef utility::pointer::shared_ptr< StoredResidueSubsetSelector const > StoredResidueSubsetSelectorCOP;

} //namespace residue_selectors
} //namespace protocols


#endif
