// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/chains_util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan, Jared Adolf-Bryfogle

#ifndef INCLUDED_core_pose_chains_util_hh
#define INCLUDED_core_pose_chains_util_hh

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

typedef std::set< int > Jumps;

/// @brief Retrieves jump information from <pose>, storing the result in <jumps>.
/// Jumps are keyed by their jump id.
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
void jumps_from_pose(core::pose::Pose const & pose, Jumps & jumps);

/// @brief get Conformation chain number -> PDBInfo chain mapping
/// @remarks Any chains whose PDBInfo chain records are marked entirely as
///  PDBInfo::empty_record() will be mapped to that character.
/// @return the mapping if PDBInfo available and chains appear consistent,
///  otherwise prints a warning and returns a default mapping (1=A, 2=B_
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
std::map< core::Size, char > conf2pdb_chain( core::pose::Pose const & pose );

/// @brief Get all the chain numbers from conformation
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// @details This is a rather silly function, as it will just return a vector
/// with entries from 1 to pose->num_chains() (as chains numbers are sequential starting from 1
utility::vector1< core::Size > get_chains( core::pose::Pose const & pose );

/// @brief compute last residue number of a chain
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
//// @details This is mostly indirection to Conformation::chain_end(), but with better error checking
core::Size chain_end_res( Pose const & pose, core::Size const chain );

/// @brief compute last residue numbers of all chains
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
//// @details This is mostly an indirection to Conformation::chain_endings(), though with better handling of the last residue
utility::vector1< core::Size > chain_end_res( Pose const & pose );

/// @brief Compute uniq chains in a complex, based on sequence identity
/// @details Returns a vector of pose length with true/false of uniq chain
///    true is unique, false is not
utility::vector1< bool > compute_unique_chains( Pose & pose );

/// @brief renumber PDBInfo based on Conformation chains; each chain starts from 1
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// @param[in,out] pose The Pose to modify.
/// @param[in] fix_chains If true, the procedure will attempt to fix any empty record
///  characters it finds in the PDBInfo. (default true)
/// @param[in] start_from_existing_numbering If true, will attempt to start each
///  chain from the existing numbering in the PDBInfo.  E.g. if the first residue
///  of chain 2 in the Conformation is 27, then the renumbering of the chain in
///  PDBInfo will start from 27. (default true)
/// @param[in] keep_insertion_codes If true, will maintain insertion codes and
///  will not increment the pdb residue numbering for those residues.  This means
///  new numbering with insertion codes will only reflect properly if the
///  old numbering included the base numbering of the insertion code residues,
///  i.e. 100 100A 100B and not just 100A 100B (with 100 never appearing).
///  (default false)
/// @param[in] rotate_chain_ids If true, allows support for more than 26 pdb chains
///  by rotating [A,Z] continuously.  WARNING: This will break the assumption
///  made by the PDBPoseMap that each pdb chain id is unique, so make sure you
///  are not using the PDBPoseMap feature downstream in your code path without
///  corrections! (default false)
/// @remarks If fixing chains and there is only one chain and the PDBInfo exists
///  but all records are marked as empty, will renumber and set the PDBInfo chain
///  to 'A'.
/// @return true if renumbering successful, false otherwise
bool renumber_pdbinfo_based_on_conf_chains(
	core::pose::Pose & pose,
	bool fix_chains = true,
	bool const start_from_existing_numbering = true,
	bool const keep_insertion_codes = false,
	bool const rotate_chain_ids = false
);

/// @brief Does the pose have a residue with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
bool
has_chain(std::string const & chain, core::pose::Pose const & pose);

/// @brief Does the pose have a residue with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
bool
has_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Does the pose have a residue with the given chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
bool
has_chain(core::Size chain_id, core::pose::Pose const & pose);

/// @brief Get all chain numbers for the residues with the given chain letters
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned chain numbers are in sorted order
utility::vector1<core::Size>
get_chain_ids_from_chains(utility::vector1<std::string> const & chains, core::pose::Pose const & pose);

/// @brief Get all chain numbers for the residues with the given chain letters
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned chain numbers are in sorted order
utility::vector1<core::Size>
get_chain_ids_from_chains(utility::vector1<char> const & chains, core::pose::Pose const & pose);

/// @brief Get all chain numbers for the residues with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned chain numbers are in sorted order
utility::vector1<core::Size>
get_chain_ids_from_chain(std::string const & chain, core::pose::Pose const & pose);

/// @brief Get all chain numbers for the residues with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned chain numbers are in sorted order
utility::vector1<core::Size>
get_chain_ids_from_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Attempt to get the chain number which correspond to the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// If the chain letter corresponds to more than one chain letter, raise an error
core::Size
get_chain_id_from_chain(std::string const & chain, core::pose::Pose const & pose);

/// @brief Attempt to get the chain number which correspond to the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// If the chain letter corresponds to more than one chain letter, raise an error
core::Size
get_chain_id_from_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Get the chain letter for the first residue in a given chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// Keep in mind that not all residues with the given chain number will necessarily have the returned chain letter
char
get_chain_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose);

/// @brief Attempt to get jump IDs which correspond to the given chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jumps here are the jumps which are *directly* upstream of a residue in the given chains,
/// (i.e. a residue on the given chain number is built directly from the jump) rather than logically upstream.
/// Return all jumps which build a given chain.
/// If no jumps directly builds the given chains (unlikely), return an empty set.
std::set<core::Size>
get_jump_ids_from_chain_ids(std::set<core::Size> const & chain_ids, core::pose::Pose const & pose);

/// @brief Attempt to get the jump number which correspond to the given chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jump here is the jump which is *directly* upstream of a residue in the given chain,
/// (i.e. a residue on the given chain number is built directly from the jump) rather than logically upstream.
/// If there's more than one jump which builds the given chain, return the smallest numbered jump
/// (even if it's not the jump which best partions the chain on the FoldTree).
/// If no jump directly builds the chain (unlikely), hard exit.
core::Size
get_jump_id_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose);

/// @brief Get all the jump numbers for the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jumps here are the jumps which are *directly* upstream of a residue with a given chain letter,
/// (i.e. a residue with the given chain letter is built directly from the jump) rather than logically upstream.
/// Return all jumps which build residues with the given chain letter.
/// If no jumps directly builds residues with the given chain letters, return an empty vector.
///
/// The returned jump numbers are in sorted order
utility::vector1<core::Size>
get_jump_ids_from_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Get all the jump numbers for the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jumps here are the jumps which are *directly* upstream of a residue with a given chain letter,
/// (i.e. a residue with the given chain letter is built directly from the jump) rather than logically upstream.
/// Return all jumps which build residues with the given chain letter.
/// If no jumps directly builds residues with the given chain letters, return an empty vector.
///
/// The returned jump numbers are in sorted order
utility::vector1<core::Size>
get_jump_ids_from_chain(std::string const & chain, core::pose::Pose const & pose);

/// @brief Get the jump number for the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jump here is the jumps which is *directly* upstream of a residue with a given chain letter,
/// (i.e. a residue with the given chain letter is built directly from the jump) rather than logically upstream.
/// If there's more than one jump which builds residues with the given chain letter, return the smallest numbered jump
/// (even if it's not the jump which best partions the chain on the FoldTree).
/// If no jump directly builds the chain, hard exit.
core::Size
get_jump_id_from_chain(std::string const & chain, core::pose::Pose const & pose);

/// @brief Get the jump number for the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The jump here is the jumps which is *directly* upstream of a residue with a given chain letter,
/// (i.e. a residue with the given chain letter is built directly from the jump) rather than logically upstream.
/// If there's more than one jump which builds residues with the given chain letter, return the smallest numbered jump
/// (even if it's not the jump which best partions the chain on the FoldTree).
/// If no jump directly builds the chain, hard exit.
core::Size
get_jump_id_from_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Get the chain ID of the residue directly built by the given jump.
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// Keep in mind that not every residue with the returned chain ID will be downstream of this jump.
core::Size
get_chain_id_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose);

/// @brief Get the chain letter of the chain built by the given jump.
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// Keep in mind that not every residue with the returned chain ID will be downstream of this jump.
char
get_chain_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose);

/// @brief Get a vector of all residues numbers which are represented by this chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned residue numbers are in sorted order
utility::vector1<core::Size>
get_resnums_for_chain( core::pose::Pose const & pose, char chain );

/// @brief Get a vector of all residues numbers which are represented by this chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
/// The returned residue numbers are in sorted order
utility::vector1<core::Size>
get_resnums_for_chain_id( core::pose::Pose const & pose, core::Size chain_id );

/// @brief Get all residues which correspond to the given chain number
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
core::conformation::ResidueCOPs
get_chain_residues(core::pose::Pose const & pose, core::Size chain_id);

/// @brief Get all residues which correspond to the given chain numbers
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
core::conformation::ResidueCOPs
get_residues_from_chains(core::pose::Pose const & pose, utility::vector1<core::Size> const & chain_ids);

/// @brief Does this residue number have this chain letter?
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
bool res_in_chain( core::pose::Pose const & pose, core::Size resnum, std::string const & chain );

/// @brief Get a value representing the position of all the atoms for residues with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
core::Size
get_hash_from_chain( char const & chain, core::pose::Pose const & pose, std::string const & extra_label="" );

/// @brief Get a value representing the position of all the atoms for residues which don't have the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
core::Size
get_hash_excluding_chain( char const & chain, core::pose::Pose const & pose, std::string const & extra_label="" );

/// @brief Get a value representing the position of all the atoms for residues with the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
std::string
get_sha1_hash_from_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label="");

/// @brief Get a value representing the position of all the atoms for residues with the given chain letters
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
std::string
get_sha1_hash_from_chains(utility::vector1< std::string > const & chains, core::pose::Pose const & pose, std::string const & extra_label="");

/// @brief Get a value representing the position of all the atoms for residues which don't have the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
std::string
get_sha1_hash_excluding_chain(char const & chain, core::pose::Pose const & pose, std::string const & extra_label="");

/// @brief Get a value representing the position of all the atoms for residues which don't have the given chain letter
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
std::string
get_sha1_hash_excluding_chains(utility::vector1< std::string > const & chains, core::pose::Pose const & pose, std::string const & extra_label="");


} // pose
} // core

#endif // INCLUDED_core_pose_chains_util_hh
