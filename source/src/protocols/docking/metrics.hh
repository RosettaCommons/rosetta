// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file metrics
/// @brief protocols that are specific to docking low resolution
/// @details
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_docking_metrics_hh
#define INCLUDED_protocols_docking_metrics_hh

// Unit Headers

// Package Headers
#include <protocols/docking/types.hh>
// Project Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers

// Numeric Headers and ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace docking {

/// @brief Calculates the difference in energy between the inputted complex, and the complex with the two partners at 500A from each other
core::Real calc_interaction_energy( const core::pose::Pose & pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps );  //@TODO These poses should be PoseCAPs!

/// @brief Calculates C-alpha RMSD of the smaller partner after superposition of the larger partner
core::Real calc_Lrmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps);

/// @brief Calculates C-alpha RMSD of the superimposed larger partner
/// @details The larger partner is the fixed partner; this is mostly for
///			membrane proteins, where we move the partners apart, then do relax
///			on both partners
core::Real calc_P1rmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps);

/// @brief Calculates C-alpha RMSD of the superimposed smaller partner
/// @details This calculates the RMSD of the smaller partner after superimposing it
///			 The computed RMSD is not from the docking, but from the flexible
///			 backbone relax
core::Real calc_P2rmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps);

/// @brief Calculates the all-atom RMSD of all residues within 5A of the interface at superposition along those same atoms
core::Real calc_Irmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps );

/// @brief calcluates the CA-atom RMSD of all residues within 5A of the interface at superposition along those same atoms
core::Real calc_CA_Irmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps );

/// @brief Calculates the fraction of native contacts recovered between the input and native pose.  A native-contact is defined
/// as defined by a residue from one partner being within 5A of a residue from the other partner in the native structure
core::Real calc_Fnat( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps );

/// @brief Calculates the fraction of non-native contacts recovered between the input and native pose.  A native-contact is defined
/// as defined by a residue from one partner being within 5A of a residue from the other partner in the native structure. Fnonnat = Nnon-native-contact/Nnative_contact
core::Real calc_Fnonnat( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps );

/// @brief Calculates the fraction of native contacts from the given native contacts pair list
core::Real calc_Fnat( const core::pose::Pose & pose, std::string const& list_file, DockJumps const movable_jumps );

/// @brief Calculates the fraction of non-native contacts from the given native contacts pari list
core::Real calc_Fnonnat( const core::pose::Pose & pose, std::string const& list_file, DockJumps const movable_jumps );

// @brief Determines if two residues are in contact within a supplied cutoff distance
bool calc_res_contact( core::conformation::ResidueOP rsd1, core::conformation::ResidueOP rsd2, core::Real dist_cutoff);

}// docking
}// protocols

#endif
