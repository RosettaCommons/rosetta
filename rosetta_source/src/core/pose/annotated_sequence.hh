// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/annotated_sequence.cc
/// @brief  utility functions for making poses from sequences
/// @author P. Douglas Renfrew
/// @author Sam Deluca

#ifndef INCLUDED_core_pose_annotated_sequence_hh
#define INCLUDED_core_pose_annotated_sequence_hh

// Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// AUTO-REMOVED #include <string>

#include <utility/vector1.hh>


namespace core {
namespace pose {

/// @brief return of list of ResidueTypes corresponding to an annotated protein sequence
/// @param[in] sequence_in an annotated sequence
/// @param[in] residue_set the desired residue set
/// @param[in] auto_termini mark position 1, last_residue with lower, upper termini; default true
chemical::ResidueTypeCAPs residue_types_from_sequence(
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini = true
);


/// @brief Creates a Pose from the annotated protein sequence  <sequence>
/// with ResidueTypeSet  <residue_set>  and stores it in  <pose>
/// @note: any existing data in  <pose>  is cleared, auto_termini
/// mark position 1, last_residue with lower, upper termini; default true
///
/// example(s):
///     make_pose_from_sequence(pose,"THANKSEVAN","fa_standard")
/// See also:
///     Pose
///     PDBInfo
///     pose_from_pdb
///     pose_from_rcsb
///     pose_from_sequence
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini = true
);


/// @brief Creates a Pose from the annotated protein sequence  <sequence>
/// with the desired  <type_set_name>  and stores it in  <pose>
/// @note: any existing data in  <pose>  is cleared, auto_termini
/// mark position 1, last_residue with lower, upper termini; default true
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence,
	std::string const & type_set_name,
	//chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini = true
);


/// @brief Returns the oneletter_sequence that corresponds to the given
/// annotated sequence.
std::string annotated_to_oneletter_sequence(
	std::string const & annotated_seq
);


} // namespace core
} // namespace pose

#endif // INCLUDED_core_pose_annotated_sequence_HH
