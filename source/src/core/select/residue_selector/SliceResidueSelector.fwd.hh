// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SliceResidueSelector.fwd.hh
/// @brief  The SliceResidueSelector allows slicing of the returned values of other residue selectors
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_SliceResidueSelector_FWD_HH
#define INCLUDED_core_select_residue_selector_SliceResidueSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace residue_selector {

class SliceResidueSelector;

typedef utility::pointer::shared_ptr< SliceResidueSelector > SliceResidueSelectorOP;
typedef utility::pointer::shared_ptr< SliceResidueSelector const > SliceResidueSelectorCOP;

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
