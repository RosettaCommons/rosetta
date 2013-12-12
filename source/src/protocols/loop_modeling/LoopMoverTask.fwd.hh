// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopMoverTask_FWD_HH
#define INCLUDED_protocols_loop_modeling_LoopMoverTask_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loop_modeling {

class LoopMoverTask;

typedef utility::pointer::owning_ptr<LoopMoverTask> LoopMoverTaskOP;
typedef utility::pointer::owning_ptr<LoopMoverTask const> LoopMoverTaskCOP;

}
}

#endif

