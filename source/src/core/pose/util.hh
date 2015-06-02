// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan

#ifndef INCLUDED_core_pose_util_hh
#define INCLUDED_core_pose_util_hh

// C/C++ headers
#include <map>
#include <set>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/SequenceMapping.fwd.hh>

// Package headers
#include <core/pose/util.tmpl.hh>
#include <core/pose/MiniPose.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/id/AtomID_Map.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#include <luabind/out_value_policy.hpp>
#endif

namespace core {
namespace pose {

typedef std::set<int> Jumps;

/// @brief Append residues of pose2 to pose1.
void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain = true
);

/// @brief Append specified residues of pose2 to pose1.
void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain = true
);

/// @brief Retrieves jump information from <pose>, storing the result in <jumps>.
/// Jumps are keyed by their jump id.
void jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps);

/// @brief Removes all virtual residues from <pose>
void remove_virtual_residues(core::pose::Pose* pose);

/// @brief Updates the rigid-body transform of the specified jump in <pose>
void swap_transform(Size jump_num, const kinematics::RT& xform, Pose* pose);

/// @brief Returns true if <residue> is positionally conserved, false otherwise
bool is_position_conserved_residue(const Pose& pose, core::Size residue);

/// @brief Create a subpose of the src pose.  PDBInfo is set as NULL.
void
create_subpose(
	Pose const & src,
	utility::vector1< Size > const & positions,
	kinematics::FoldTree const & f,
	Pose & pose
);

/// @brief Create a subpose of the src pose -- figures out a reasonable fold tree.
void
pdbslice( pose::Pose & new_pose,
					pose::Pose const & pose,
					utility::vector1< Size > const & slice_res );

/// @brief Create a subpose of the src pose -- figures out a reasonable fold tree.
void
pdbslice( pose::Pose & pose,
					utility::vector1< Size > const & slice_res );

// for partition_by_jump: both new poses start residue numbering from 1 and don't keep the original numbering!
void
partition_pose_by_jump(
	pose::Pose const & src,
	int const jump_number,
	pose::Pose & partner1, // partner upstream in foldtree
	pose::Pose & partner2  // partner downstream in foldtree
);

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
// 	pose::Pose &pose
// );

utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose );

/// getters/setters for things in the Pose DataCache

// ARBITRARY_FLOAT_DATA
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const name,
	core::Real & value
);

Real getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const name );

bool
hasPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const name );

void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const name,
	core::Real value
);

void clearPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name
);

void clearPoseExtraScores(
	core::pose::Pose & pose
);

// ARBITRARY_STRING_DATA
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const name,
	std::string & value
);

void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string name,
	std::string const value
);

/// @brief Adds a VRT res to the end of the pose at the center of mass.
/// Reroots the structure on this res.
void addVirtualResAsRoot(core::pose::Pose & pose);

/// @brief Adds a virtual residue to the end of the pose at the specified location.
/// Roots the structure on this residue.
void addVirtualResAsRoot(const numeric::xyzVector<core::Real>& xyz, core::pose::Pose& pose);

/// @brief Get center of mass of a pose.
numeric::xyzVector< core::Real >
get_center_of_mass( core::pose::Pose const & pose );

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

/// @brief get Conformation chain -> PDBInfo chain mapping
/// @remarks Any chains whose PDBInfo chain records are marked entirely as
///  PDBInfo::empty_record() will be mapped to that character.  Note that
///  Conformation -> PDBInfo is always unique, but the reverse may not be true.
/// @return the mapping if PDBInfo available and chains appear consistent,
///  otherwise returns an empty mapping
std::map< int, char > conf2pdb_chain( core::pose::Pose const & pose );


/// @brief renumber PDBInfo based on Conformation chains; each chain starts from 1
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

/// @brief Returns true if the  <pose>  geometry is ideal
/// @param[in] pose The Pose to check.
/// @return true if all pose positions have ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_pose(
									 core::pose::Pose const & pose
);

/// @brief Returns true if the  <pose> geometry is ideal in position  <seqpos>
/// @param[in] pose The Pose to check.
/// @return true if position seqpos has ideal bond lengths and angles
///  up to some very small epsilon
bool is_ideal_position(
											 Size seqpos,
											 core::pose::Pose const & pose
);

/// @brief this function removes all residues from the pose which are not protein residues.  This removal includes, but is not limited to, metals, DNA, RNA, and ligands.  It will NOT remove ligands which are canonical residues (for example, if a protein binds an alanine monomer, the monomer will be untouched).
void remove_nonprotein_residues( core::pose::Pose & pose );

/// @brief this function removes all residues with both UPPER and LOWER terminus types.  This is intended for removing ligands that are canonical residues.
void remove_ligand_canonical_residues( core::pose::Pose & pose );

/// @brief this function compares pose atom coordinates for equality; it is not the == operator because it does not compare all pose data.
/// @author Steven Lewis smlewi@gmail.com
/// @param[in] lhs one pose to compare
/// @param[in] rhs one pose to compare
/// @param[in] n_dec_places number of decimal places to compare for the coordinates (remember == doesn't work for float); defaults to 3 which is PDB accuracy
bool compare_atom_coordinates(
	Pose const & lhs,
	Pose const & rhs,
	Size const n_dec_places = 3);

/// @brief this function compares poses for equality up to the
///information stored in the binary protein silent struct format.
bool compare_binary_protein_silent_struct(
	Pose const & lhs,
	Pose const & rhs);


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

id::NamedAtomID
atom_id_to_named_atom_id(
	core::id::AtomID const & atom_id,
	Pose const & pose
);

id::AtomID
named_atom_id_to_atom_id(
	core::id::NamedAtomID const & named_atom_id,
	Pose const & pose,
	bool raise_exception = true
);

id::NamedStubID
stub_id_to_named_stub_id(
	core::id::StubID const & stub_id,
	Pose const & pose
);

id::StubID
named_stub_id_to_stub_id(
	core::id::NamedStubID const & named_stub_id,
	Pose const & pose
);

///////////////////////////////////////////////////////////////////
std::string tag_from_pose( core::pose::Pose const & pose );
void tag_into_pose( core::pose::Pose & pose, std::string const & tag );

// criterion for sorting.
bool sort_pose_by_score( core::pose::PoseOP const & pose1, core::pose::PoseOP const & pose2 );

core::Real energy_from_pose(
	core::pose::Pose const & pose, core::scoring::ScoreType const & sc_type
);

core::Real energy_from_pose(
	core::pose::Pose const & pose, std::string const & sc_type
);

core::Real total_energy_from_pose( core::pose::Pose const & pose );

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose, core::Size ir, core::Size jr );

void
transfer_phi_psi( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose );

void
transfer_jumps( const core::pose::Pose& srcpose, core::pose::Pose& tgtpose);

void
replace_pose_residue_copying_existing_coordinates(
	pose::Pose & pose,
	Size const seqpos,
	core::chemical::ResidueType const & new_rsd_type
	);

/// @brief Remove variant from an existing residue.
conformation::ResidueOP remove_variant_type_from_residue(
		conformation::Residue const & old_rsd,
		core::chemical::VariantType const variant_type,
		pose::Pose const & pose );

/// @brief Construct a variant of an existing residue.
conformation::ResidueOP add_variant_type_to_residue(
		conformation::Residue const & old_rsd,
		core::chemical::VariantType const variant_type,
		pose::Pose const & pose );

/// @brief Construct a variant of an existing pose residue.
void add_variant_type_to_pose_residue(
		pose::Pose & pose,
		chemical::VariantType const variant_type,
		Size const seqpos );

/// @brief Construct a non-variant of an existing pose residue.
void remove_variant_type_from_pose_residue(
		pose::Pose & pose,
		chemical::VariantType const variant_type,
		Size const seqpos );


void
add_lower_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	);


void
add_upper_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	);


void
remove_lower_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	);


void
remove_upper_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
	);

/// @brief set up a map to look up TORSION_ID by DOF_ID (Map[DOF_ID] = TORISION_ID)
void
setup_dof_to_torsion_map(
	pose::Pose const & pose,
	id::DOF_ID_Map< id::TorsionID > & dof_map
	);


/// @brief convert from allow-bb/allow-chi MoveMap to simple DOF_ID boolean mask needed by the minimizer
void
setup_dof_mask_from_move_map(
		core::kinematics::MoveMap const & mm,
		pose::Pose const & pose,
		id::DOF_ID_Mask & dof_mask
);

bool
has_chain(std::string const & chain, core::pose::Pose const & pose);

bool
has_chain(char const & chain, core::pose::Pose const & pose);

bool
has_chain(core::Size chain_id, core::pose::Pose const & pose);

std::set<core::Size>
get_jump_ids_from_chain_ids(std::set<core::Size> const chain_ids, core::pose::Pose const & pose);

core::Size
get_jump_id_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose);

//std::set<core::Size> get_chain_ids(core::pose::Pose const & pose);

core::Size
get_chain_id_from_chain(std::string const & chain, core::pose::Pose const & pose);

core::Size
get_chain_id_from_chain(char const & chain, core::pose::Pose const & pose);

utility::vector1<core::Size>
get_chain_ids_from_chain(std::string const & chain, core::pose::Pose const & pose);

utility::vector1<core::Size>
get_chain_ids_from_chain(char const & chain, core::pose::Pose const & pose);

utility::vector1<core::Size>
get_chain_ids_from_chains(utility::vector1<std::string> const & chains, core::pose::Pose const & pose);

char
get_chain_from_chain_id(core::Size const & chain_id, core::pose::Pose const & pose);

core::Size
get_jump_id_from_chain(std::string const & chain, core::pose::Pose const & pose);

core::Size
get_jump_id_from_chain(char const & chain, core::pose::Pose const & pose);

utility::vector1<core::Size>
get_jump_ids_from_chain(char const & chain, core::pose::Pose const & pose);

utility::vector1<core::Size>
get_jump_ids_from_chain(std::string const & chain, core::pose::Pose const & pose);

core::Size
get_chain_id_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose);

char
get_chain_from_jump_id(core::Size const & jump_id, core::pose::Pose const & pose);

core::conformation::ResidueCOPs
get_chain_residues(core::pose::Pose const & pose, core::Size chain_id);

/// @brief Is residue number in this chain?
bool res_in_chain( core::pose::Pose const & pose, core::Size resnum, std::string chain );

core::Size num_heavy_atoms(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
);

core::Size num_atoms(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
);

core::Size num_hbond_acceptors(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
);

core::Size num_hbond_donors(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
);
core::Size
num_chi_angles(
	core::Size begin,
	core::Size const end,
	core::pose::Pose const & pose
);

core::Real
mass(
		core::Size begin,
		core::Size const end,
		core::pose::Pose const & pose
);

core::Size
get_hash_from_chain(char const & chain, core::pose::Pose const & pose);

core::Size
get_hash_excluding_chain(char const & chain, core::pose::Pose const & pose);

std::string get_sha1_hash_excluding_chain(char const & chain, core::pose::Pose const & pose);

/// @brief Initialize a DOF_ID_Map for a given Pose using the DOF_ID_Map's current default fill values
template< typename T >
void
initialize_dof_id_map( id::DOF_ID_Map< T > & dof_map, Pose const & pose );

/// @brief Initialize a DOF_ID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_dof_id_map( id::DOF_ID_Map< T > & dof_map, Pose const & pose, T const & value );

/// @brief returns a Distance
core::Real
pose_max_nbr_radius( pose::Pose const & pose );

/// @brief Initialize an AtomID_Map for a given Pose using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, pose::Pose const & pose );

/// @brief Initialize an AtomID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, pose::Pose const & pose, T const & value );

/// @brief Initialize an AtomID_Map for a given Conformation using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation );

/// @brief Initialize an AtomID_Map for a given Conformation using a specified fill value
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation, T const & value );

/// @brief Initialize an AtomID_Map for a given Pose using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, pose::Pose const & pose );

/// @brief Initialize an AtomID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, pose::Pose const & pose, T const & value );

/// @brief Initialize an AtomID_Map for a given Conformation using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation );

/// @brief Initialize an AtomID_Map for a given Conformation using a specified fill value
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation, T const & value );

/// @brief detect and fix disulfide bonds
void
initialize_disulfide_bonds(
	Pose & pose);

/// @brief Returns a string giving the pose's tag if there is such a thing or "UnknownTag" otherwise.
std::string extract_tag_from_pose( core::pose::Pose &pose );

/// @brief Create a sequence map of first pose onto the second, matching the PDBInfo
///    If the PDBInfo of either Pose is missing or invalid, do a simple sequence alignment matching.
core::id::SequenceMapping sequence_map_from_pdbinfo( Pose const & first, Pose const & second );
#ifdef USELUA
void lregister_util( lua_State * lstate );
#endif

/// @brief count the number of canonical residues in the pose
core::Size canonical_residue_count(core::pose::Pose const & pose);

/// @brief count the number of non-canonical residues in the pose
core::Size noncanonical_residue_count(core::pose::Pose const & pose);

/// @brief count the number of canonical amino acid atoms in the pose
core::Size canonical_atom_count(core::pose::Pose const & pose);

/// @brief count the number of non-canonical amino acids in thepose
core::Size noncanonical_atom_count(core::pose::Pose const & pose);

/// @brief count the number of non-canonical chi angles in the pose
core::Size noncanonical_chi_count(core::pose::Pose const & pose);

/// @brief Number of protein residues in the pose
/// @details No virtuals, membrane residues or embedding residues counted
core::Size nres_protein( core::pose::Pose const & pose );

numeric::xyzVector< core::Real>
center_of_mass(
	core::pose::Pose const & pose,
	int const start,
	int const stop
);

int
residue_center_of_mass(
	core::pose::Pose const & pose,
	int const start,
	int const stop
);

int
return_nearest_residue(
	core::pose::Pose const & pose,
	int const begin,
	int const end,
	core::Vector center
);

id::AtomID_Map< id::AtomID >
convert_from_std_map( std::map< id::AtomID, id::AtomID > const & atom_map, core::pose::Pose const & pose );

/// @brief Add cutpoint variants to all residues annotated as cutpoints in the pose.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose );

void
correctly_add_cutpoint_variants( core::pose::Pose & pose,
																 Size const cutpoint_res,
																 bool const check_fold_tree = true );

/// @brief returns true if the given residue in the pose is a chain ending or has upper/lower terminal variants
bool
pose_residue_is_terminal( Pose const & pose, Size const resid );

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_lower_terminus( pose::Pose const & pose, Size const resid );

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_upper_terminus( pose::Pose const & pose, Size const resid );

} // pose
} // core

#endif // INCLUDED_core_pose_util_HH
