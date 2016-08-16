// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/NeighborhoodResidueSelector.fwd.hh
/// @brief  Forward declaration of a class that selects residues from the neighborhood of a set of focus residues
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_NeighborhoodResidueSelector_FWD_HH
#define INCLUDED_core_select_residue_selector_NeighborhoodResidueSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace select {
namespace residue_selector {

class NeighborhoodResidueSelector;

typedef utility::pointer::shared_ptr< NeighborhoodResidueSelector > NeighborhoodResidueSelectorOP;
typedef utility::pointer::shared_ptr< NeighborhoodResidueSelector const > NeighborhoodResidueSelectorCOP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
