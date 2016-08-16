// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LigandMotifSearch.fwd.hh
/// @brief Forward declaration for interaction motif search protocol
/// @author mdsmith (mdwsmith@u.washington.edu)

#ifndef INCLUDED_protocols_motifs_LigandMotifSearch_FWD_HH
#define INCLUDED_protocols_motifs_LigandMotifSearch_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace motifs {

class LigandMotifSearch;
typedef utility::pointer::shared_ptr< LigandMotifSearch > LigandMotifSearchOP;
typedef utility::pointer::shared_ptr< LigandMotifSearch const > LigandMotifSearchCOP;

} // motifs
} // protocols

#endif // INCLUDED_protocols_motifs_LigandMotifSearch_fwd
