// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/Loop.fwd.hh
/// @brief  Loop class forward declarations header
/// @author Chu Wang


#ifndef INCLUDED_protocols_loops_Loop_FWD_HH
#define INCLUDED_protocols_loops_Loop_FWD_HH

// Unit header
#include <utility/pointer/owning_ptr.hh>

// Project header
#include <core/types.hh>

// Utility Header
#include <utility/vector1.hh>

namespace protocols {
namespace loops {

class Loop;
typedef utility::pointer::shared_ptr< Loop > LoopOP;
typedef utility::pointer::shared_ptr< Loop const > LoopCOP;

struct SerializedLoop;

typedef utility::vector1< SerializedLoop > SerializedLoopList;

//typedef utility::pointer::shared_ptr< SerializedLoopList > SerializedLoopListOP;
//typedef utility::pointer::shared_ptr< SerializedLoopList const > SerializedLoopListCOP;

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_Loop_FWD_HH
