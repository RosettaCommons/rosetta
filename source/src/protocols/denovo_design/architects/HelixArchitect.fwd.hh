// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/HelixArchitect.fwd.hh
/// @brief Architect for helices
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_denovo_design_architects_HelixArchitect_fwd_hh
#define INCLUDED_protocols_denovo_design_architects_HelixArchitect_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace denovo_design {
namespace architects {

class HelixArchitect;

typedef utility::pointer::shared_ptr< HelixArchitect > HelixArchitectOP;
typedef utility::pointer::shared_ptr< HelixArchitect const > HelixArchitectCOP;



} //protocols
} //denovo_design
} //architects


#endif //INCLUDED_protocols_denovo_design_architects_HelixArchitect_fwd_hh





