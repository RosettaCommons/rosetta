// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/GlycanInfoMover.fwd.hh
/// @brief Simple class for outputting glycan information. Currently, it simply prints the information.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_GlycanInfoMover_fwd_hh
#define INCLUDED_protocols_analysis_GlycanInfoMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace analysis {

class GlycanInfoMover;

typedef utility::pointer::shared_ptr< GlycanInfoMover > GlycanInfoMoverOP;
typedef utility::pointer::shared_ptr< GlycanInfoMover const > GlycanInfoMoverCOP;

} //protocols
} //analysis

#endif //INCLUDED_protocols_analysis_GlycanInfoMover_fwd_hh
