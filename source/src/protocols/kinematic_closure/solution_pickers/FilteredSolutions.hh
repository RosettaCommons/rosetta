// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_solution_pickers_FilteredSolutions_HH
#define INCLUDED_protocols_kinematic_closure_solution_pickers_FilteredSolutions_HH

// Unit headers
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace solution_pickers {

/// @brief Apply rama and bump checks to quickly filter out bad solutions.
///
/// @details The rama check only considers the pivot torsions, since it is 
/// assumed that the non-pivot torsions were set as desired.  This check  
/// employs a Monte Carlo acceptance step, which means that the pivot torsions 
/// will really be sampled from a rama distribution (ignoring the inherent 
/// geometric biases in the closure algorithm) when this check is enabled.  
/// However, note that the rama distribution will be double-counted if the 
/// score function also contains a rama term, which it usually does.
///
/// The bump check checks for clashes between the N, CA, C, O, and CB atoms of 
/// every residue in the loop versus every other residue in the protein.  This 
/// filter is O(n^2) and is much slower than the constant time rama check.  The 
/// rama check is also very selective; it usually filters out more than 90% of 
/// proposed solutions.  For these reasons, it is important that the rama check 
/// be run before the bump check.  You can disable the rama check, but this 
/// would probably lead to a noticeable drop in performance.

class FilteredSolutions : public SolutionPicker {

public:
	/// @brief Constructor which can enable or disable any of the filters used by 
	/// this algorithm.
	FilteredSolutions(
			bool check_rama=true,
			bool check_overlap=true,
			bool be_lenient=false);

public:
	/// @copydoc SolutionPicker::pick_and_apply
	bool pick_and_apply(Pose & pose, SolutionList const & solutions);

public:
	/// @brief Enable the rama check.
	void check_rama() { check_rama_ = true; }

	/// @brief Enable the bump check.
	void check_overlap() { check_overlap_ = true; }

	/// @brief Make the bump check more lenient.
	void be_lenient() { be_lenient_ = true; }

	/// @brief Disable the rama check.
	void dont_check_rama() { check_rama_ = false; }

	/// @brief Disable the bump check.
	void dont_check_overlap() { check_overlap_ = false; }

	/// @brief Make the bump check more lenient.
	void dont_be_lenient() { be_lenient_ = false; }

private:
	bool check_rama_, check_overlap_, be_lenient_;

};

}
}
}

#endif

