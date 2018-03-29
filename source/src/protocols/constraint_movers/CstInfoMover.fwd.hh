// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_movers/CstInfoMover.fwd.hh
/// @brief A Mover to output information about constraints
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_protocols_constraint_movers_CstInfoMover_fwd_hh
#define INCLUDED_protocols_constraint_movers_CstInfoMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace constraint_movers {

class CstInfoMover;

typedef utility::pointer::shared_ptr< CstInfoMover > CstInfoMoverOP;
typedef utility::pointer::shared_ptr< CstInfoMover const > CstInfoMoverCOP;


} //protocols
} //constraint_movers


#endif //INCLUDED_protocols_constraint_movers_CstInfoMover_fwd_hh





