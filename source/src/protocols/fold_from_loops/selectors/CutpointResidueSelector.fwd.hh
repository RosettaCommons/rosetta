// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/CutpointResidueSelector.fwd.hh
/// @brief  Selects those residues that are cutpoints of the FoldTree
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_CutpointResidueSelector_FWD_HH
#define INCLUDED_protocols_fold_from_loops_CutpointResidueSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

// Project Headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace fold_from_loops {
namespace selectors {

class CutpointResidueSelector;

typedef utility::pointer::shared_ptr< CutpointResidueSelector > CutpointResidueSelectorOP;
typedef utility::pointer::shared_ptr< CutpointResidueSelector const > CutpointResidueSelectorCOP;

} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#endif
