// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/LigandConformer.fwd.hh
/// @brief  Forward declaration for a class to hold the geometry of a ligand conformation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_LigandConformer_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_LigandConformer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class LigandConformer;

typedef utility::pointer::shared_ptr< LigandConformer > LigandConformerOP;
typedef utility::pointer::shared_ptr< LigandConformer const > LigandConformerCOP;

}
}
}

#endif
