// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/minimization_packing/symmetry/SymMinMover.fwd.hh
/// @brief  MinMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_fwd_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace minimization_packing {
namespace symmetry {

//Forwards and OP typedefs
class SymMinMover;
typedef utility::pointer::shared_ptr< SymMinMover > SymMinMoverOP;
typedef utility::pointer::shared_ptr< SymMinMover const > SymMinMoverCOP;

} // symmetry
}//moves
}//protocols

#endif //INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_FWD_HH
