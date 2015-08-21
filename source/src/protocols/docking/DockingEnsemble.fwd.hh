// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/DockingEnsemble.fwd.hh
/// @brief  DockingEnsemble forward declarations header
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_docking_DockingEnsemble_FWD_HH
#define INCLUDED_protocols_docking_DockingEnsemble_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace docking {

//Forwards and OP typedefs
class DockingEnsemble;
typedef utility::pointer::shared_ptr< DockingEnsemble > DockingEnsembleOP;
typedef utility::pointer::shared_ptr< DockingEnsemble const > DockingEnsembleCOP;

}//docking
}//protocols

#endif //INCLUDED_protocols_moves_DockingEnsemble_FWD_HH
