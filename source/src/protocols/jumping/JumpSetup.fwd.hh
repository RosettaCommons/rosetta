// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jumping_JumpSetup_fwd_hh
#define INCLUDED_protocols_jumping_JumpSetup_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace jumping {

// Forward
class BaseJumpSetup;
class JumpSetup;
class JumpSelector;

// Types
typedef  utility::pointer::shared_ptr< BaseJumpSetup >  BaseJumpSetupOP;
typedef  utility::pointer::shared_ptr< BaseJumpSetup const >  BaseJumpSetupCOP;

typedef  utility::pointer::shared_ptr< JumpSetup >  JumpSetupOP;
typedef  utility::pointer::shared_ptr< JumpSetup const >  JumpSetupCOP;

typedef  utility::pointer::shared_ptr< JumpSelector >  JumpSelectorOP;
typedef  utility::pointer::shared_ptr< JumpSelector const >  JumpSelectorCOP;


} // namespace kinematics
} // namespace core

#endif
