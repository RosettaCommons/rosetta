// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/calc_taskop_movers/ForceDisulfidesMover.fwd.hh
/// @brief  ForceDisulfidesMover forward declarations header
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_moves_ForceDisulfidesMover_fwd_hh
#define INCLUDED_protocols_simple_moves_ForceDisulfidesMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace calc_taskop_movers {

//Forwards and OP typedefs
class ForceDisulfidesMover;
typedef utility::pointer::shared_ptr< ForceDisulfidesMover > ForceDisulfideMoverOP;
typedef utility::pointer::shared_ptr< ForceDisulfidesMover const > ForceDisulfideMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_ForceDisulfidesMover_FWD_HH
