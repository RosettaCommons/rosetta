// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan, Jared Adolf-Bryfogle

#ifndef INCLUDED_core_pose_util_hh
#define INCLUDED_core_pose_util_hh

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
#include <tuple>
#include <set>


namespace core {
namespace pose {

/// @brief Updates the rigid-body transform of the specified jump in <pose>
void swap_transform(Size jump_num, kinematics::RT const & xform, Pose & pose);

/// @brief Returns true if <residue> is positionally conserved, false otherwise
/// Based on the POSITION_CONSERVED_RESIDUES annotation stored in the Pose DataCache
bool is_position_conserved_residue(Pose const & pose, core::Size residue);

void set_reasonable_fold_tree( core::pose::Pose & pose );


/// @brief Return the appropritate ResidueType for the virtual residue for the i
/// "mode" (fullatom, centroid ...) the pose is in.
core::chemical::ResidueTypeCOP
virtual_type_for_pose(core::pose::Pose const & pose);

/// @brief Adds a VRT res to the end of the pose at the center of mass.
/// Reroots the structure on this res.
void addVirtualResAsRoot(core::pose::Pose & pose);

/// @brief Adds a virtual residue to the end of the pose at the specified location.
/// Roots the structure on this residue.
void addVirtualResAsRoot(const numeric::xyzVector<core::Real>& xyz, core::pose::Pose& pose);

/// @brief Removes all virtual residues from <pose>
void remove_virtual_residues(core::pose::Pose & pose);

/// @brief Get center of mass of a pose.
///
/// This computes an equally-weighted, all-(non-virtual)-heavy atom center.
numeric::xyzVector< core::Real >
get_center_of_mass( core::pose::Pose const & pose );

/// @brief Repair pdbinfo of inserted residues that may have blank chain and zero
/// seqpos. Assumes insertions only occur _after_ a residue.
void fix_pdbinfo_damaged_by_insertion(
	core::pose::Pose & pose
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

/// @brief Return the residue type in the correct
/// "mode" (fullatom, centroid ...) the pose is in.
core::chemical::ResidueTypeCOP
get_restype_for_pose(core::pose::Pose const & pose, std::string const & name);

/// @brief Return the residue type in the passed mode,
/// respecting any modification that pose may make.
core::chemical::ResidueTypeCOP
get_restype_for_pose(core::pose::Pose const & pose, std::string const & name, core::chemical::TypeSetMode mode);

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
	Pose & pose
);

/// @brief detect and fix disulfide bonds
void
initialize_disulfide_bonds(
	Pose & pose,
	io::StructFileRep const & fd
);

/// @brief Create a sequence map of first pose onto the second, matching the PDBInfo
///    If the PDBInfo of either Pose is missing or invalid, do a simple sequence alignment matching.
core::id::SequenceMapping sequence_map_from_pdbinfo( Pose const & first, Pose const & second );

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

/// @brief Get the center of the indicated residues
///
/// WARNING: Despite the name, this function only calculates with a single coordinate per residue
/// (the Calpha/neighbor atom)
numeric::xyzVector< core::Real>
center_of_mass(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & residues
);

/// @brief Get the center of the indicated residues
///
/// WARNING: Despite the name, this function only calculates with a single coordinate per residue
/// (the Calpha/neighbor atom)
numeric::xyzVector< core::Real>
center_of_mass(
	core::pose::Pose const & pose,
	int const start,
	int const stop
);

/// @brief Get the center of the indicated residues
///
/// This computes an equally-weighted, all-atom (including virtuals and hydrogens) center
core::Vector
all_atom_center(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residues
);

int
residue_center_of_mass(
	pose::Pose const & pose,
	utility::vector1< bool > residues
);

int
return_nearest_residue(
	pose::Pose const & pose,
	utility::vector1< bool > const & residues,
	Vector center
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

/// @brief Create std::map from PDBPoseMap
std::map< std::string, core::Size > get_pdb2pose_numbering_as_stdmap ( core::pose::Pose const & pose );

/// @brief Create a chemical bond from lower to upper residue across CUTPOINT_LOWER/CUTPOINT_UPPER.
/// @details This will prevent steric repulsion.
/// @param[in] pose The pose to modify.
/// @param[in] cutpoint_res The index of the CUTPOINT_LOWER residue.
/// @param[in] next_res_in The index of the CUTPOINT_UPPER residue.  If not provided, or if set to 0, this defaults
/// to the cutpoint_res + 1 residue.  Must be specified for cyclic geometry.
void
declare_cutpoint_chemical_bond( core::pose::Pose & pose, Size const cutpoint_res, Size const next_res_in = 0 );

/// @brief Given a pose and a position that may or may not be CUTPOINT_UPPER or CUTPOINT_LOWER, determine whether this
/// position has either of these variant types, and if it does, determine whether it's connected to anything.  If it is,
/// update the C-OVL1-OVL2 bond lengths and bond angle (for CUTPOINT_LOWER) or OVU1-N bond length (for CUTPOINT_UPPER) to
/// match any potentially non-ideal geometry in the residue to which it's bonded.
/// @details Requires a little bit of special-casing for gamma-amino acids.  Throws an exception if the residue to which
/// a CUTPOINT_LOWER is bonded does not have an "N" and a "CA" or "C4".  Safe to call repeatedly, or if cutpoint variant
/// types are absent; in these cases, the function does nothing.
/// @note By default, this function calls itself again once on residues to which this residue is connected, to update their
/// geometry.  Set recurse=false to disable this.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
update_cutpoint_virtual_atoms_if_connected( core::pose::Pose & pose, core::Size const cutpoint_res, bool recurse = true );

void
get_constraints_from_link_records( core::pose::Pose & pose, io::StructFileRep const & sfr );

/// @brief Convert PDB numbering to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pdb_to_pose( pose::Pose const & pose, utility::vector1< int > const & pdb_res );

/// @brief Convert PDB numbering/chain to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pdb_to_pose( pose::Pose const & pose, std::tuple< utility::vector1< int >, utility::vector1<char>, utility::vector1<std::string> > const & pdb_res );

/// @brief Convert PDB numbering to pose numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
Size pdb_to_pose( pose::Pose const & pose, int const res_num, char const chain = ' ' );

/// @brief Convert pose numbering to pdb numbering. Must exist somewhere else, but I couldn't find it. -- rhiju
utility::vector1< Size > pose_to_pdb( pose::Pose const & pose, utility::vector1< Size > const & pose_res );

/// @brief  Is the query atom in this pose residue axial or equatorial to the given ring or neither?
chemical::rings::AxEqDesignation is_atom_axial_or_equatorial_to_ring(
	Pose const & pose,
	uint seqpos,
	uint query_atom,
	utility::vector1< uint > const & ring_atoms );

/// @brief  Is the query atom in this pose axial or equatorial to the given ring or neither?
/*chemical::rings::AxEqDesignation is_atom_axial_or_equatorial_to_ring(
Pose const & pose,
id::AtomID const & query_atom,
utility::vector1< id::AtomID > const & ring_atoms );*/


/// @brief  Is the query atom in this pose residue axial or equatorial or neither?
//chemical::rings::AxEqDesignation is_atom_axial_or_equatorial( Pose const & pose, uint seqpos, uint query_atom );

/// @brief  Is the query atom in this pose axial or equatorial or neither?
//chemical::rings::AxEqDesignation is_atom_axial_or_equatorial( Pose const & pose, id::AtomID const & query_atom );

/// @brief Set bfactors in a pose PDBInfo
void
set_bfactors_from_atom_id_map(Pose & pose, id::AtomID_Map< Real > const & bfactors);



///Set the BB torsion, phi, psi, omega (see core::types).
/// Works with carbohydrates.
/// Think about moving this to pose itself.
void
set_bb_torsion( uint torsion_id, Pose & pose, core::Size sequence_position, core::Angle new_angle);

///@brief Get a particular backbone torsion, phi, psi, omega (see core::types)
/// Works with carbohydrates.
/// Think about moving this to pose itself.
core::Angle
get_bb_torsion( uint torsion_id, Pose const & pose, core::Size sequence_position );


// Stepwise

bool
just_modeling_RNA( std::string const & sequence );

bool
stepwise_addable_pose_residue( core::Size const n /*in pose numbering*/,
	core::pose::Pose const & pose );

bool
stepwise_addable_residue( core::Size const n /*in full model numbering*/,
	std::map< core::Size, std::string > const & non_standard_residue_map );


////////////////////////////////////////////////////////////////////////////////////////////////
bool
effective_lower_terminus_based_on_working_res( Size const i,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & cutpoint_open_in_full_model );

bool
effective_upper_terminus_based_on_working_res( Size const i,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	Size const nres_full);

bool
definite_terminal_root( utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & res_list,
	Size const nres,
	Size const i );

bool
definite_terminal_root( pose::Pose const & pose, Size const i );

Size
get_definite_terminal_root( pose::Pose const & pose,
	utility::vector1< Size > const & partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res );

Size
get_definite_terminal_root( pose::Pose const & pose,
	utility::vector1< Size > const & partition_res /* should not be empty */ );

utility::vector1< Size >
reorder_root_partition_res(
	utility::vector1< Size > const & root_partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */ );

void
reroot( pose::Pose & pose,
	utility::vector1< Size > const & root_partition_res /* should not be empty */,
	utility::vector1< Size > const & res_list,
	utility::vector1< Size > const & preferred_root_res /* can be empty */,
	utility::vector1< Size > const & fixed_domain_map /* 0 in free; 1,2,... for separate fixed domains */,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & working_res );

} // pose
} // core

#endif // INCLUDED_core_pose_util_HH
