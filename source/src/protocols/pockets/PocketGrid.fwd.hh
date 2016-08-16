// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/PocketGrid.fwd.hh
/// @brief  protocols::pockets::PocketGrid forward declarations header
/// @author David Johnson
/// @author Ragul Gowthaman


#ifndef INCLUDED_protocols_pockets_PocketGrid_fwd_hh
#define INCLUDED_protocols_pockets_PocketGrid_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace pockets {

// Forward
class PCluster;
class PClusterSet;
class PocketGrid;
class EggshellGrid;

// Types
typedef utility::pointer::shared_ptr< PocketGrid >  PocketGridOP;
typedef utility::pointer::shared_ptr< PocketGrid const >  PocketGridCOP;


} // namespace pockets
} // namespace protocols

#endif
