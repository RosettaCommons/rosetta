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
/// @author Labonte (carbohydrate versions)

#ifndef INCLUDED_core_pose_annotated_sequence_hh
#define INCLUDED_core_pose_annotated_sequence_hh

// Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AA.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace pose {

/// @brief Parse the input annotated sequence
void parse_sequence(
	std::string const & sequence_in,
	utility::vector1< std::string > & fullname_list,
	std::vector< Size > & oneletter_to_fullname_index,
	std::string & one_letter_sequence
);

/// @brief Get the real length of a annotated sequence
Size get_sequence_len( std::string const & sequence_in );


/// @brief return a list of ResidueTypes corresponding to an annotated protein sequence
/// @param[in] sequence_in an annotated sequence
/// @param[in] residue_set the desired residue set
/// @param[in] auto_termini mark position 1, last_residue with lower, upper termini; default true
chemical::ResidueTypeCOPs residue_types_from_sequence(
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini = true
);


/// @brief Return a list of carbohydrate ResidueTypes corresponding to an annotated, linear, IUPAC polysaccharide
/// sequence.
chemical::ResidueTypeCOPs residue_types_from_saccharide_sequence( std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set );


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
	chemical::ResidueTypeCOPs requested_types,
	bool const auto_termini = true
);

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
	bool const auto_termini = true
);


/// @brief Create a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with ResidueTypeSet
/// <residue_set> and store it in <pose>.
void make_pose_from_saccharide_sequence( pose::Pose & pose,
	std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini = true );

/// @brief Create a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
/// <type_set_name> and store it in <pose>.
void make_pose_from_saccharide_sequence( pose::Pose & pose,
	std::string const & sequence,
	std::string const & type_set_name = "fa_standard",
	bool const auto_termini = true );

/// @brief Return a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
/// <type_set_name>.
pose::PoseOP pose_from_saccharide_sequence( std::string const & sequence,
	std::string const & type_set_name = "fa_standard",
	bool const auto_termini = true );


/// @brief Returns the oneletter_sequence that corresponds to the given
/// annotated sequence.
std::string annotated_to_oneletter_sequence(
	std::string const & annotated_seq
);

/// @brief use efficient residue type finder to find simplest residue type with this AA & requested termini.
core::chemical::ResidueTypeCOP
get_rsd_type_from_aa( chemical::ResidueTypeSet const & residue_set,
	core::chemical::AA const & my_aa, bool const & is_lower_terminus, bool const & is_upper_terminus );

core::chemical::ResidueTypeCOP
get_rsd_type_from_aa_legacy( chemical::ResidueTypeSet const & residue_set,
	core::chemical::AA const & my_aa, bool const & is_lower_terminus, bool const & is_upper_terminus );

core::chemical::ResidueTypeCOP
handle_disulfide_variant_name_for_backwards_compatibility( core::chemical::ResidueTypeSet const & residue_set,
	std::string const &fullname );


} // namespace core
} // namespace pose

#endif // INCLUDED_core_pose_annotated_sequence_HH
