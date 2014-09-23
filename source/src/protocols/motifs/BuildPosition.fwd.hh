// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BuildPosition.fwd.hh
/// @brief Forward declaration for position where motifs are built in MotifSearch
/// @author sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_BuildPosition_fwd_hh
#define INCLUDED_protocols_motifs_BuildPosition_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace motifs {

class BuildPosition;
typedef utility::pointer::shared_ptr< BuildPosition > BuildPositionOP;
typedef utility::pointer::shared_ptr< BuildPosition const > BuildPositionCOP;
typedef utility::vector1< BuildPositionOP > BuildPositionOPs;
typedef utility::vector1< BuildPositionCOP > BuildPositionCOPs;

} // motifs
} // protocols

#endif // INCLUDED_protocols_motifs_BuildPosition_fwd
