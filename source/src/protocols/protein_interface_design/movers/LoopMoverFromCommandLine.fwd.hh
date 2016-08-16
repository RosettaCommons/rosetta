// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/ProteinInterfaceDesign/movers/LoopRemodel.fwd.hh
/// @brief  LoopRemodel forward declarations header

#ifndef INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLine_FWD_HH
#define INCLUDED_protocols_protein_interface_design_movers_LoopMoverFromCommandLine_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

//Forwards and OP typedefs
class LoopMoverFromCommandLine;
typedef utility::pointer::shared_ptr< LoopMoverFromCommandLine > LoopMoverFromCommandLineOP;
typedef utility::pointer::shared_ptr< LoopMoverFromCommandLine const > LoopMoverFromCommandLineCOP;

}//moves
}//protocols
}

#endif
