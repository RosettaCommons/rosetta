// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ProteinResidueSelector.fwd.hh
/// @brief  Selects those residues that are amino acids
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_ProteinResidueSelector_FWD_HH
#define INCLUDED_protocols_fold_from_loops_ProteinResidueSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace fold_from_loops {
namespace selectors {

class ProteinResidueSelector;

typedef utility::pointer::shared_ptr< ProteinResidueSelector > ProteinResidueSelectorOP;
typedef utility::pointer::shared_ptr< ProteinResidueSelector const > ProteinResidueSelectorCOP;

} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#endif
