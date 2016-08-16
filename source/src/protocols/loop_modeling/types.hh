// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_types_HH
#define INCLUDED_protocols_loop_modeling_types_HH

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loop_modeling/LoopMover.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace loop_modeling {

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::LoopsOP;
using protocols::loops::LoopsCOP;
using utility::vector1;

typedef vector1<Size> IndexList;
typedef vector1<LoopMoverOP> LoopMovers;

}
}

#endif

