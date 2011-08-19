// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MoverStub.fwd.hh
/// @brief  MoverStub forward declarations header
/// @author

#ifndef INCLUDED_protocols_moves_MoverStub_FWD_HH
#define INCLUDED_protocols_moves_MoverStub_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace moves{

//Forwards and OP typedefs
class MoverStub;
typedef utility::pointer::owning_ptr< MoverStub > MoverStubOP;
typedef utility::pointer::owning_ptr< MoverStub const > MoverStubCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_MoverStub_FWD_HH
