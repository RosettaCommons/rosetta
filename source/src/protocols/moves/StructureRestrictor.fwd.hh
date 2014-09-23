// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/StructureRestrictor.fwd.hh
/// @brief  lookup relevant chains for a structure in a table.
/// @author Matthew O'Meara

/// This should probably be a pilot app, but the way Rosetta Scripts
/// is set up, it can't be in the pilot apps

#ifndef INCLUDED_protocols_moves_StructureRestrictor_fwd_hh
#define INCLUDED_protocols_moves_StructureRestrictor_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace moves{

class StructureRestrictor;
typedef utility::pointer::shared_ptr< StructureRestrictor > StructureRestrictorOP;
typedef utility::pointer::shared_ptr< StructureRestrictor const > StructureRestrictorCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_StructureRestrictor_FWD_HH
