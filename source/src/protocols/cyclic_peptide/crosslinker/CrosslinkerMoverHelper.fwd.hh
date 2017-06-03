// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.fwd.hh
/// @brief A base class for helper objects that the CrosslinkerMover uses to set up specific types
/// of threefold linkers.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_CrosslinkerMoverHelper_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_CrosslinkerMoverHelper_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

class CrosslinkerMoverHelper;

typedef utility::pointer::shared_ptr< CrosslinkerMoverHelper > CrosslinkerMoverHelperOP;
typedef utility::pointer::shared_ptr< CrosslinkerMoverHelper const > CrosslinkerMoverHelperCOP;

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_crosslinker_CrosslinkerMoverHelper_fwd_hh
