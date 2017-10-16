// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/extra_pose_info_util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan, Jared Adolf-Bryfogle

#ifndef INCLUDED_core_pose_extra_pose_info_util_hh
#define INCLUDED_core_pose_extra_pose_info_util_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/MiniPose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rings/AxEqDesignation.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <map>
#include <set>


namespace core {
namespace pose {

/// @brief Analyzes  <pose>  residue phi/psi sets and guesses the secondary
/// structure, ideally dssp should be used for that
void
set_ss_from_phipsi(
	pose::Pose & pose
);

// /// @brief Analyses the pose in terms of phi/psi and guesses at the secondary
// /// structure - ideally dssp should be used for that
// void
// set_ss_from_phipsi_dssp(
//  pose::Pose &pose
// );

utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose, std::string const & filename );
utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose );

/// getters/setters for things in the Pose DataCache

/// @brief return bool is T/F for whether the requested datum exists.  "value" is the data, pass-by-ref.
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	core::Real & value
);

/// @brief return value is ExtraScore if exist, runtime_assert if it doesn't exist
Real getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name );

/// @brief does this ExtraScore exist?
bool
hasPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name );

void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	core::Real value
);

void clearPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name
);

void clearPoseExtraScores(
	core::pose::Pose & pose
);

/// @brief return bool is T/F for whether the requested datum exists.  "value" is the data, pass-by-ref.
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	std::string & value
);

void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	std::string const & value
);

/// @brief Adds a key-value pair to the STRING_MAP in the Pose DataCache. If
/// there is no STRING_MAP in the DataCache, one is created.
void add_comment(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
);

/// @brief Attempts to access the entry in the STRING_MAP associated with the
/// given key. If an entry for the key exists, the value associated with the key
/// is put into val, and this function returns true. Otherwise, this function
/// returns false and val left unmodified.
bool get_comment(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
);

/// @brief Deletes the entry in the STRING_MAP associated with the
/// given key.
void delete_comment(
	core::pose::Pose & pose,
	std::string const & key
);

/// @brief Gets a map< string, string > representing comments about the Pose in
/// the form of key-value pairs.
std::map< std::string, std::string > get_all_comments(
	core::pose::Pose const & pose
);

/// @brief Dumps a pdb with comments at end of file


/// @brief Sets a PDB-style REMARK entry in the Pose.
/// @details This is different from a comment in its interpretation by the
/// silent-file output machinery. A REMARK is written on its own separate line
/// in the output silent-file, while a comment is written as part of the Pose
/// SCORE: lines.
void add_score_line_string(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
);

bool get_score_line_string(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
);

/// @brief Gets a map< string, string > representing score_line_strings about the Pose in
/// the form of key-value pairs.
std::map< std::string, std::string > get_all_score_line_strings(
	core::pose::Pose const & pose
);

/// @brief  Reads the comments from the pdb file and adds it into comments
void read_comment_pdb(
	std::string const &file_name,
	core::pose::Pose  & pose
);
/// @brief  dumps pose+ comments to pdb file
void dump_comment_pdb(
	std::string const &file_name,
	core::pose::Pose const& pose
);

std::string tag_from_pose( core::pose::Pose const & pose );
/// @brief Returns a string giving the pose's tag if there is such a thing or "UnknownTag" otherwise.
std::string extract_tag_from_pose( core::pose::Pose &pose );
void tag_into_pose( core::pose::Pose & pose, std::string const & tag );

} // pose
} // core

#endif // INCLUDED_core_pose_extra_pose_info_util_HH
