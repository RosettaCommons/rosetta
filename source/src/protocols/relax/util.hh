// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/relax/util.hh
/// @brief initialization protocols for relax and utility functions
///
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_util_hh
#define INCLUDED_protocols_relax_util_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/relax/RelaxProtocolBase.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace relax {

/// @brief (pba)  read in membrane topology
void setup_membrane_topology( core::pose::Pose & pose, std::string spanfile );

/// @brief Make DNA Rigid in movemap if found in pose
void make_dna_rigid(core::pose::Pose & pose, core::kinematics::MoveMap & mm);

/// @brief Set energy method options for DNA-DNA.  Any other optimization for DNA relax should go here.
void setup_for_dna(core::scoring::ScoreFunction & scorefxn);

void relax_pose( core::pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag );

RelaxProtocolBaseOP generate_relax_from_cmd( bool NULL_if_no_cmd = false );


/// @brief Rebuild all polymer bond-dependent atoms in the pose.
/// @details Rewritten from an earlier, uglier function that manually placed H atom in non-pro.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).
void
fixH (
	core::pose::Pose & pose
);

/// @brief Add terminal peptide bond and cyclic constraints to pose.
/// @author Gaurav Bhardwaj, Baker laboratory
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu) -- moved from protocols/forge/methods/util.hh/cc
void
cyclize_pose(core::pose::Pose & pose);

}
} // protocols

#endif
