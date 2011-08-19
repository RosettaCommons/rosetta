// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange



#ifndef INCLUDED_protocols_evaluation_util_hh
#define INCLUDED_protocols_evaluation_util_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>

// Package Headers
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/rms_util.hh>

// ObjexxFCL Headers
//// C++ headers
#include <list>

//Auto Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>


namespace protocols {
namespace evaluation {

///@brief register cmd-line options relevant for evaluators...
void register_options();

void read_common_evaluator_options( MetaPoseEvaluator& );

//@detail find residues that don't have missing density
void find_existing_residues(  core::pose::PoseCOP pose, std::string tag, core::scoring::ResidueSelection& selection );
void invert_include_residues( core::Size nres, core::scoring::ResidueSelectionVector const& include_list, core::scoring::ResidueSelectionVector& exclude_list );
void evaluate_pose( core::pose::Pose& pose, PoseEvaluator& eval, std::ostream& );

// this function will return a bunch of "loops" that refer to residues that are considered part of the core:
// not scored are loops with 4 or more residues, short helices (<=5) that terminate a loop are not scored, too
void define_scorable_core_from_secondary_structure( core::fragment::SecondaryStructure const&, protocols::loops::Loops& score_core );

}
}
#endif
