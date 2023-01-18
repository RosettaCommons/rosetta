// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/PyRotamerEliminator.fwd.hh
/// @brief  Wrappers that allow the user to eliminate rotamers from pyrosetta
/// @author Jack Maguire


#ifndef INCLUDED_core_pack_rotamer_set_PyRotamerEliminator_fwd_hh
#define INCLUDED_core_pack_rotamer_set_PyRotamerEliminator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {

class PyRotamerEliminator;
using PyRotamerEliminatorOP = utility::pointer::shared_ptr< PyRotamerEliminator >;
using PyRotamerEliminatorCOP = utility::pointer::shared_ptr< PyRotamerEliminator const >;

class PyRotamerEliminatorTaskOperation;
using PyRotamerEliminatorTaskOperationOP = utility::pointer::shared_ptr< PyRotamerEliminatorTaskOperation >;
using PyRotamerEliminatorTaskOperationCOP = utility::pointer::shared_ptr< PyRotamerEliminatorTaskOperation const >;

}
}
}

#endif
