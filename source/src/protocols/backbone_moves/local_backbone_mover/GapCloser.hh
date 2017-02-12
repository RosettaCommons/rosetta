// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/GapCloser.hh
/// @brief GapCloser closes the gaps after moving the free peptide.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_hh

// Project headers
#include <protocols/backbone_moves/local_backbone_mover/GapCloser.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/types.hh>
#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/GapSolutionPicker.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

///@brief GapCloser closes the gaps after moving the free peptide.
class GapCloser : public utility::pointer::ReferenceCount {
	
	friend class ::LocalBackboneMoverTests;

public:

	GapCloser();
	GapCloser(GapCloser const & src);

	virtual ~GapCloser();

	GapCloserOP
	clone() const;

	/// @brief Set the solution picker
	void set_solution_picker(gap_solution_pickers::GapSolutionPickerOP solution_picker){
		solution_picker_ = solution_picker;
	}
	
	/// @brief Find solutions for the two gaps
	void solve_gaps(FreePeptide &free_peptide);

	/// @brief Return true if the gap problem is solved
	bool gap_solved(){return gap_solved_;}

	/// @brief Apply closure to a pose
	/// @detail Exception will be raised if the closure problems haven't been solved
	void apply_closure(core::pose::Pose &pose, FreePeptide &free_peptide);

private: // Member functions

	/// @brief Solve a gap closing problem
	bool solve_a_gap(FreePeptide &free_peptide, Size pivot, vector1<vector1<Real> > &pivot_torsions);

	/// @brief Pick a pair of solutions
	void pick_solutions(Size &index1, Size &index2, core::pose::Pose &pose, FreePeptide &free_peptide);

private: // Member data

	bool gap_solved_ = false;

	vector1<vector1<Real> > pivot_torsions1_; // Solution of torsions for the first gap
	vector1<vector1<Real> > pivot_torsions2_; // Solution of torsions for the second gap

	gap_solution_pickers::GapSolutionPickerOP solution_picker_ = nullptr;
};


} //protocols
} //backbone_moves
} //local_backbone_mover



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_GapCloser_hh





