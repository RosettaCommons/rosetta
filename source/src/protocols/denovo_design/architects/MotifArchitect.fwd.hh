// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/MotifArchitect.fwd.hh
/// @brief Designs topologies
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_MotifArchitect_fwd_hh
#define INCLUDED_protocols_denovo_design_architects_MotifArchitect_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// Forward
namespace protocols {
namespace denovo_design {
namespace architects {

class MotifArchitect;

typedef utility::pointer::shared_ptr< MotifArchitect > MotifArchitectOP;
typedef utility::pointer::shared_ptr< MotifArchitect const > MotifArchitectCOP;

typedef utility::vector1< MotifArchitectOP > MotifArchitectOPs;
typedef utility::vector1< MotifArchitectCOP > MotifArchitectCOPs;

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_MotifArchitect_fwd_hh

