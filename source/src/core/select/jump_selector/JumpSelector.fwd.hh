// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelector.fwd.hh
/// @brief  Forward declaration of a class that identifies a subset of residues from a Pose
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_JumpSelector_FWD_HH
#define INCLUDED_core_select_jump_selector_JumpSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace select {
namespace jump_selector {

typedef utility::vector1< bool > JumpSubset;
typedef utility::pointer::shared_ptr< JumpSubset > JumpSubsetOP;
typedef utility::pointer::shared_ptr< JumpSubset const > JumpSubsetCOP;
typedef utility::pointer::weak_ptr< JumpSubset > JumpSubsetAP;
typedef utility::pointer::weak_ptr< JumpSubset const > JumpSubsetCAP;

class JumpSelector;

typedef utility::pointer::shared_ptr< JumpSelector > JumpSelectorOP;
typedef utility::pointer::shared_ptr< JumpSelector const > JumpSelectorCOP;
typedef utility::pointer::weak_ptr< JumpSelector > JumpSelectorAP;
typedef utility::pointer::weak_ptr< JumpSelector const > JumpSelectorCAP;

} //namespace jump_selector
} //namespace select
} //namespace core


#endif
