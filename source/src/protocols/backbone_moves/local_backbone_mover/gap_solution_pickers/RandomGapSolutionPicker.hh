// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/RandomGapSolutionPicker.hh
/// @brief Pick a solution randomly.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_RandomGapSolutionPicker_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_RandomGapSolutionPicker_hh

#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/GapSolutionPicker.hh>
#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/RandomGapSolutionPicker.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace gap_solution_pickers {

///@brief Pick a solution randomly.
class RandomGapSolutionPicker : public GapSolutionPicker {

public:

	RandomGapSolutionPicker();

	virtual ~RandomGapSolutionPicker();

	virtual Size pick(core::pose::Pose &pose, FreePeptide &free_peptide,
		vector1<vector1<Real> > pivot_torsions, Size pivot);

private:

};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //gap_solution_pickers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_RandomGapSolutionPicker_hh





