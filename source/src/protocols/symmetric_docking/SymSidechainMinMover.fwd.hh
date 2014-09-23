// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RestrictTaskForDocking.fwd.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_symmetric_docking_SymSidechainMinMover_fwd_hh
#define INCLUDED_protocols_symmetric_docking_SymSidechainMinMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace symmetric_docking {

class SymSidechainMinMover;
typedef utility::pointer::shared_ptr< SymSidechainMinMover > SymSidechainMinMoverOP;
typedef utility::pointer::shared_ptr< SymSidechainMinMover > SymSidechainMinMoverCOP;

class SymInterfaceSidechainMinMover;
typedef utility::pointer::shared_ptr< SymInterfaceSidechainMinMover > SymInterfaceSidechainMinMoverOP;
typedef utility::pointer::shared_ptr< SymInterfaceSidechainMinMover > SymInterfaceSidechainMinMoverCOP;

}
}

#endif
