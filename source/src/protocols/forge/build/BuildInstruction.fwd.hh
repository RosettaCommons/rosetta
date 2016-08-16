// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/BuildInstruction.fwd.hh
/// @brief forward declaration for BuildInstruction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_BuildInstruction_fwd_hh
#define INCLUDED_protocols_forge_build_BuildInstruction_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief forward declaration for BuildInstruction
class BuildInstruction;


/// @brief owning pointer for BuildInstruction
typedef utility::pointer::shared_ptr< BuildInstruction > BuildInstructionOP;


/// @brief owning pointer for const BuildInstruction
typedef utility::pointer::shared_ptr< BuildInstruction const > BuildInstructionCOP;


/// @brief access pointer for BuildInstruction
typedef utility::pointer::weak_ptr< BuildInstruction > BuildInstructionAP;


/// @brief access pointer for const BuildInstruction
typedef utility::pointer::weak_ptr< BuildInstruction const > BuildInstructionCAP;


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_BuildInstruction_FWD_HH */
