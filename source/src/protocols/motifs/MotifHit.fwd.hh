// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifHit.fwd.hh
/// @brief Forward declaration for class that holds information about a motif in the context of the search
/// @author sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_MotifHit_fwd_hh
#define INCLUDED_protocols_motifs_MotifHit_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace motifs {

class MotifHit;
typedef utility::pointer::shared_ptr< MotifHit > MotifHitOP;
typedef utility::pointer::shared_ptr< MotifHit const > MotifHitCOP;
typedef utility::vector1< MotifHitOP > MotifHitOPs;
typedef utility::vector1< MotifHitCOP > MotifHitCOPs;

} // motifs
} // protocols

#endif // INCLUDED_protocols_motifs_MotifHit_fwd
