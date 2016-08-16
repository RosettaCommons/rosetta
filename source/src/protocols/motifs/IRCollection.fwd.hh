// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IRCollection.fwd.hh
/// @brief forward declaration for sets of interaction motifis between residues
/// @author havranek


#ifndef INCLUDED_protocols_motifs_IRCollection_fwd_hh
#define INCLUDED_protocols_motifs_IRCollection_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace motifs {

class IRCollection;
typedef utility::pointer::shared_ptr< IRCollection > IRCollectionOP;
typedef utility::pointer::shared_ptr< IRCollection const > IRCollectionCOP;
typedef utility::vector1< IRCollectionOP > IRCollectionOPs;

}
}

#endif // INCLUDED_protocols_motifs_IRCollection_fwd
