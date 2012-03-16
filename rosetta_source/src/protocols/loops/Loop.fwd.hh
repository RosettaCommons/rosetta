// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/Loop.fwd.hh
/// @brief  Loop class forward declarations header
/// @author Chu Wang


#ifndef INCLUDED_protocols_loops_Loop_FWD_HH
#define INCLUDED_protocols_loops_Loop_FWD_HH

// Unit header
#include <utility/pointer/owning_ptr.hh>

// Project header
#include <core/types.hh>

namespace protocols {
namespace loops {

class Loop;
typedef utility::pointer::owning_ptr< Loop > LoopOP;
typedef utility::pointer::owning_ptr< Loop const > LoopCOP;

struct SerializedLoop {
	core::Size start;
	core::Size stop;
	core::Size cut;
	core::Real skip_rate;
	bool extended;
};

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_Loop_FWD_HH
