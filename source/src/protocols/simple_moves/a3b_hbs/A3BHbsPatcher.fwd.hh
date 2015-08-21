// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/a3b_hbs/A3BHbsPatcher.fwd.hh
/// @brief  A3BHbsPatcher forward declarations header
/// @author Andy Watkins, amw579@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsPatcher_fwd_hh
#define INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsPatcher_fwd_hh
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace simple_moves {
namespace a3b_hbs {


//Forwards and OP typedefs
class A3BHbsPatcher;
typedef utility::pointer::shared_ptr< A3BHbsPatcher > A3BHbsPatcherOP;
typedef utility::pointer::shared_ptr< A3BHbsPatcher const > A3BHbsPatcherCOP;

}//hbs
}//simple_moves
}//protocols

#endif //INCLUDED_protocols_simple_moves_hbs_HbsPatcher_fwd_hh
