// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Chu Wang


#ifndef INCLUDED_protocols_jumping_ResiduePairJumpSetup_fwd_hh
#define INCLUDED_protocols_jumping_ResiduePairJumpSetup_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace jumping {

// Forward
class ResiduePairJumpSetupSingle;

// Types
typedef  utility::pointer::shared_ptr< ResiduePairJumpSetupSingle >  ResiduePairJumpSetupSingleOP;
typedef  utility::pointer::shared_ptr< ResiduePairJumpSetupSingle const >  ResiduePairJumpSetupSingleCOP;

// Forward
class ResiduePairJumpSetup;

// Types
typedef  utility::pointer::shared_ptr< ResiduePairJumpSetup >  ResiduePairJumpSetupOP;
typedef  utility::pointer::shared_ptr< ResiduePairJumpSetup const >  ResiduePairJumpSetupCOP;

} // namespace kinematics
} // namespace core

#endif
