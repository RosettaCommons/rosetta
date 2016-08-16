// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SecondaryStructureSelector.fwd.hh
/// @brief  Forward declaration of a class to select residues based on their secondary structure
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_select_residue_selector_SecondaryStructureSelector_FWD_HH
#define INCLUDED_core_select_residue_selector_SecondaryStructureSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

// Project Headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace core {
namespace select {
namespace residue_selector {

class SecondaryStructureSelector;

typedef utility::pointer::shared_ptr< SecondaryStructureSelector > SecondaryStructureSelectorOP;
typedef utility::pointer::shared_ptr< SecondaryStructureSelector const > SecondaryStructureSelectorCOP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
