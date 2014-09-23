// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/simple_moves/symmetry/SymMinMover.fwd.hh
/// @brief  MinMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_simple_moves_symmetry_SymMinMover_fwd_hh
#define INCLUDED_protocols_simple_moves_symmetry_SymMinMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace simple_moves{
namespace symmetry {

//Forwards and OP typedefs
class SymMinMover;
typedef utility::pointer::shared_ptr< SymMinMover > SymMinMoverOP;
typedef utility::pointer::shared_ptr< SymMinMover const > SymMinMoverCOP;

} // symmetry
}//moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_symmetry_SymMinMover_FWD_HH
