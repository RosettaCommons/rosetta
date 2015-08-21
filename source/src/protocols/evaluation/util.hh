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
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_util_hh
#define INCLUDED_protocols_evaluation_util_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/rms_util.hh>

// ObjexxFCL Headers
//// C++ headers
#include <list>

#include <core/fragment/SecondaryStructure.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

/// @brief register cmd-line options relevant for evaluators...
void register_options();

//@detail find residues that don't have missing density
void find_existing_residues(  core::pose::PoseCOP pose, std::string tag, core::scoring::ResidueSelection& selection );
void invert_include_residues( core::Size nres, core::scoring::ResidueSelectionVector const& include_list, core::scoring::ResidueSelectionVector& exclude_list );
void evaluate_pose( core::pose::Pose& pose, PoseEvaluator& eval, std::ostream& );


}
}
#endif
