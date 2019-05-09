// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/AqueousPoreFinder.fwd.hh
/// @brief Compute the center, major axis and minor axis of an ellipsoidal aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_AqueousPoreFinder_fwd_hh
#define INCLUDED_protocols_membrane_AqueousPoreFinder_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace membrane {

class AqueousPoreFinder;

typedef utility::pointer::shared_ptr< AqueousPoreFinder > AqueousPoreFinderOP;
typedef utility::pointer::shared_ptr< AqueousPoreFinder const > AqueousPoreFinderCOP;

} //protocols
} //membrane

#endif //INCLUDED_protocols_membrane_AqueousPoreFinder_fwd_hh
