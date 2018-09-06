// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file   core/pack/task/residue_selector/NativeSelector.fwd.hh
/// @brief  A ResidueSelector that applies a given residue selector to the native pose.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_residue_selectors_NativeSelector_fwd_hh
#define INCLUDED_protocols_residue_selectors_NativeSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace residue_selectors {

class NativeSelector;

typedef utility::pointer::shared_ptr< NativeSelector > NativeSelectorOP;
typedef utility::pointer::shared_ptr< NativeSelector const > NativeSelectorCOP;

} //protocols
} //residue_selectors


#endif //INCLUDED_protocols_residue_selectors_NativeSelector_fwd_hh





