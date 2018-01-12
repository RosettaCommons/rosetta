// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/minimization_packing/MinMover.fwd.hh
/// @brief  MinMover forward declarations header
/// @author

#ifndef INCLUDED_protocols_minimization_packing_SaneMinMover_fwd_hh
#define INCLUDED_protocols_minimization_packing_SaneMinMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace minimization_packing {

//Forwards and OP typedefs
class SaneMinMover;
typedef utility::pointer::shared_ptr< SaneMinMover > SaneMinMoverOP;
typedef utility::pointer::shared_ptr< SaneMinMover const > SaneMinMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_minimization_packing_SaneMinMover_FWD_HH
