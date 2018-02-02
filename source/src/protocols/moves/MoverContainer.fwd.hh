// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverContainer.fwd.hh
/// @brief  MoverContainer forward declarations header
/// @author

#ifndef INCLUDED_protocols_moves_MoverContainer_fwd_hh
#define INCLUDED_protocols_moves_MoverContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

//Forwards and OP typedefs
class MoverContainer;
typedef utility::pointer::shared_ptr< MoverContainer > MoverContainerOP;
typedef utility::pointer::shared_ptr< MoverContainer const > MoverContainerCOP;

class SequenceMover;
typedef utility::pointer::shared_ptr< SequenceMover > SequenceMoverOP;
typedef utility::pointer::shared_ptr< SequenceMover const > SequenceMoverCOP;

class RandomMover;
typedef utility::pointer::shared_ptr< RandomMover > RandomMoverOP;
typedef utility::pointer::shared_ptr< RandomMover const > RandomMoverCOP;

class CycleMover;
typedef utility::pointer::shared_ptr< CycleMover > CycleMoverOP;
typedef utility::pointer::shared_ptr< CycleMover const > CycleMoverCOP;

class SwitchMover;
typedef utility::pointer::shared_ptr< SwitchMover > SwitchMoverOP;
typedef utility::pointer::shared_ptr< SwitchMover const > SwitchMoverCOP;

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_MoverContainer_FWD_HH
