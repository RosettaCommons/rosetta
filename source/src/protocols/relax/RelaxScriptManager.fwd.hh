// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/RelaxScriptManager.fwd.hh
/// @brief A singleton class for managing relax scripts, to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_relax_RelaxScriptManager_fwd_hh
#define INCLUDED_protocols_relax_RelaxScriptManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace relax {

class EnergyMapContainer;
class RelaxScriptManager;
class RelaxScriptFileContents;

typedef utility::pointer::shared_ptr< EnergyMapContainer > EnergyMapContainerOP;
typedef utility::pointer::shared_ptr< EnergyMapContainer const > EnergyMapContainerCOP;
typedef utility::pointer::shared_ptr< RelaxScriptManager > RelaxScriptManagerOP;
typedef utility::pointer::shared_ptr< RelaxScriptManager const > RelaxScriptManagerCOP;
typedef utility::pointer::shared_ptr< RelaxScriptFileContents > RelaxScriptFileContentsOP;
typedef utility::pointer::shared_ptr< RelaxScriptFileContents const > RelaxScriptFileContentsCOP;

} //protocols
} //relax

#endif //INCLUDED_protocols_relax_RelaxScriptManager_fwd_hh
