// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/GlycanSequonsSelector.fwd.hh
/// @brief Find glycosylation sequons in pose
/// @author raemisch (raemisch@scripps.edu)

#ifndef INCLUDED_core_select_residue_selector_GlycanSequonsSelector_fwd_hh
#define INCLUDED_core_select_residue_selector_GlycanSequonsSelector_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace select {
namespace residue_selector {

class GlycanSequonsSelector;

typedef utility::pointer::shared_ptr< GlycanSequonsSelector > GlycanSequonsSelectorOP;
typedef utility::pointer::shared_ptr< GlycanSequonsSelector const > GlycanSequonsSelectorCOP;

} //core
} //select
} //residue_selector

#endif //INCLUDED_core_select_residue_selector_GlycanSequonsSelector_fwd_hh
