// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/annotated_sequence.cc
/// @brief  utility functions for making poses from sequences
/// @author P. Douglas Renfrew
/// @author Sam Deluca
/// @author Labonte (carbohydrate versions)


#ifndef INCLUDED_core_pose_annotated_sequence_hh
#define INCLUDED_core_pose_annotated_sequence_hh

// Package header
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AA.hh>

// Utility header
#include <utility/vector1.hh>

// C++ header
#include <string>


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


/// @brief Return a list of carbohydrate ResidueTypes corresponding to an annotated, linear, IUPAC polysaccharide sequence.
/// @param[in] <sequence>: an annotated IUPAC polysaccharide sequence,
/// e.g., "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp"
/// @param[in] <residue_set>: the desired residue set
/// @return    a 1-indexed vector of ResidueType owning pointers, in Rosetta order.
/// That is, the first residue will be the reducing end (the last residue in the sequence),
/// and all branches will be consecutive, with earlier branches earlier in set of residue types.
/// @details   Format for <sequence>:\n
/// Prefixes apply to the residue to which they are attached, below indicated by residue n.\n
/// Residues are listed from N to 1, where N is the total number of residues in the saccharide.\n
/// The sequence is parsed by reading to the next hyphen, so hyphens are crucial.\n
/// Linkage indication: "(a->x)-" specifies the linkage of residue n, where a is the anomeric carbon number of residue
/// (n+1) and x is the oxygen number of residue n.  The first residue listed in the annotated sequence (residue N)
/// need not have the linkage prefix.  A ->4) ResidueType will automatically be assigned by default if not specified.\n
/// Anomer indication: The strings "alpha-" or "beta-" are supplied next, which determines the stereochemistry of the
/// anomeric carbon of the residue to which it is prefixed.  An alpha ResidueType will automatically be assigned by
/// default.\n
/// Stereochemical indication: "L-" or "D-" specifies whether residue n is an L- or D-sugar.  The default is "D-".\n
/// 3-Letter code: A three letter code (in sentence case) MUST be supplied next.  This specifies the "base sugar name",
/// e.g., Glc is for glucose.  (A list of all recognized 3-letter codes for sugars can be found in the database.)\n
/// 1-Letter suffix: If no suffix follows, residue n will be linear.  If a letter is present, it indicates the ring
/// size, where "f" is furanose, "p" is pyranose, and "s" is septanose.\n
/// Branches are indicated using nested brackets and are best explained by example:\n
/// beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc is:\n
/// beta-D-Galp-(1->4)-D-GlcpNAc\n
///                       |\n
///     alpha-L-Fucp-(1->3)
/// @note make_pose_from_saccharide_sequence() will generate a pose with a proper lower terminus.
/// glycosylate_pose() will append the fragment by bond.
chemical::ResidueTypeCOPs residue_types_from_saccharide_sequence( std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set );

/// @brief  Append an empty or current Pose with saccharide residues, building branches as necessary.
void append_pose_with_glycan_residues(
	pose::Pose & pose,
	chemical::ResidueTypeCOPs residue_types,
	core::uint resnum_to_be_appended = 0 );


/// @brief Creates a Pose from the annotated protein sequence  <sequence>
/// with ResidueTypeSet  <residue_set>  and stores it in  <pose>
/// @note: any existing data in  <pose>  is cleared, auto_termini
/// mark position 1, last_residue with lower, upper termini; default true
///
/// example(s):
///     make_pose_from_sequence(pose,"THANKSEVAN",core::chemical::FA_STANDARD)
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

void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence,
	chemical::ResidueTypeSetCOP residue_set,
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
	bool const auto_termini = true,
	bool const idealize_linkages = true );

/// @brief Create a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
/// <type_set_name> and store it in <pose>.
void make_pose_from_saccharide_sequence( pose::Pose & pose,
	std::string const & sequence,
	std::string const & type_set_name = "fa_standard",
	bool const auto_termini = true,
	bool const idealize_linkages = true );

/// @brief Return a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
/// <type_set_name>.
pose::PoseOP pose_from_saccharide_sequence( std::string const & sequence,
	std::string const & type_set_name = "fa_standard",
	bool const auto_termini = true,
	bool const idealize_linkages = true );


/// @brief Returns the oneletter_sequence that corresponds to the given
/// annotated sequence.
std::string annotated_to_oneletter_sequence( std::string const & annotated_seq );

/// @brief use efficient residue type finder to find simplest residue type with this AA & requested termini.
core::chemical::ResidueTypeCOP
get_rsd_type_from_aa( chemical::ResidueTypeSet const & residue_set,
	core::chemical::AA const & my_aa, bool const & is_lower_terminus, bool const & is_upper_terminus );

} // namespace core
} // namespace pose

#endif // INCLUDED_core_pose_annotated_sequence_HH
