// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ncbb/a3b_hbs/A3BHbsMover.fwd.hh
/// @brief  HbsMover forward declarations header
/// @author Kevin Drew, kdrew@nyu.edu and Andy Watkins, amw579@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsMover_fwd_hh
#define INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsMover_fwd_hh
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace simple_moves{
namespace a3b_hbs{


//Forwards and OP typedefs
class A3BHbsMover;
typedef utility::pointer::shared_ptr< A3BHbsMover > A3BHbsMoverOP;
typedef utility::pointer::shared_ptr< A3BHbsMover const > A3BHbsMoverCOP;

}//hbs
}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_hbs_HbsMover_fwd_hh
