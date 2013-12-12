// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_solution_pickers_RandomSolutions_HH
#define INCLUDED_protocols_kinematic_closure_solution_pickers_RandomSolutions_HH

// Unit headers
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.hh>
#include <protocols/kinematic_closure/solution_pickers/RandomSolutions.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace solution_pickers {

/// @brief Pick a completely random solution to return.
class RandomSolutions : public SolutionPicker {

public:
	/// @copydoc SolutionPicker::pick_and_apply
	bool pick_and_apply(Pose & pose, SolutionList const & solutions);

};

}
}
}

#endif

