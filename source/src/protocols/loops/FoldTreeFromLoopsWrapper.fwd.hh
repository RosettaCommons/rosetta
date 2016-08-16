// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/protein_interface_design/movers/FoldTreeFromLoops.fwd.hh
/// @brief  FoldTreeFromLoops forward declarations header
/// @author Sarel

#ifndef INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_FWD_HH
#define INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {

//Forwards and OP typedefs
class FoldTreeFromLoops;
typedef utility::pointer::shared_ptr< FoldTreeFromLoops > FoldTreeFromLoopsOP;
typedef utility::pointer::shared_ptr< FoldTreeFromLoops const > FoldTreeFromLoopsCOP;

}//loops
}//protocols

#endif //INCLUDED_protocols_loops_FoldTreeFromLoops_FWD_HH
