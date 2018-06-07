// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/util.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_util_HH
#define INCLUDED_core_conformation_util_HH


// Unit headers
#include <core/conformation/Conformation.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/rings/AxEqDesignation.hh>

#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>

#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/RT.fwd.hh>

// C++ headers
#include <iosfwd>

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#include <string>
#endif

namespace core {
namespace conformation {

void
orient_residue_for_ideal_bond(
	Residue & moving_rsd,
	chemical::ResidueConnection const & moving_connection,
	Residue const & fixed_rsd,
	chemical::ResidueConnection const & fixed_connection,
	Conformation const & conformation,
	bool lookup_bond_length = false
);

/// @brief  Sets the two bond angles and the bond length across the junction, rebuilds dependent atoms (eg O,H)
void
insert_ideal_bonds_at_polymer_junction(
	Size const seqpos,
	Conformation & conformation
);

void
insert_ideal_mainchain_bonds(
	Size const seqpos,
	Conformation & conformation
);


/// @brief  Idealize backbone and sidechain at seqpos
void
idealize_position(
	Size const seqpos,
	Conformation & conformation
);


/// @brief  Return true if position contains an ideal geometry up to some epsilon
///
/// @params seqpos - sequence position
/// @params conformation - conformation object
/// @params theta_epsilon - permitted deviation from ideal bond angles, in Radians
/// @params D_epsilon - permitted deviation from ideal bond length
///
/// @remarks conformation is needed for context of polymer nbrs
bool
is_ideal_position( // Barak 6/30/09
	Size const seqpos,
	Conformation const& conformation,
	Real theta_epsilon = 0.005, // ~0.29 degrees
	Real D_epsilon = 0.02
);


/// @brief  Fills coords of target_rsd with coords from source_rsd of same atom_name, rebuilds others.
/// @details  If preserve_only_sidechain_dihedrals is true, then this function only copies mainchain coordinates,
/// and rebuilds all sidechain coordinates from scratch, setting side-chain dihedrals based on the source residue.
/// Otherwise, if false, it copies all the atoms that it can from the source residue, then rebuilds the rest.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
copy_residue_coordinates_and_rebuild_missing_atoms(
	Residue const & source_rsd,
	Residue & target_rsd,
	Conformation const & conformation,
	bool const preserve_only_sidechain_dihedrals
);


/// @brief  Fills coords of target_rsd with coords from source_rsd of same atom_name, rebuilds others.
///
void
copy_residue_coordinates_and_rebuild_missing_atoms(
	Residue const & source_rsd,
	Residue & target_rsd,
	Conformation const & conformation
);

/// @brief Given two residues, check that they are compatible types to be connected via a cutpoint.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
check_good_cutpoint_neighbour(
	core::conformation::Residue const &thisres,
	core::conformation::Residue const &other_res
);

/// @brief Given a conformation and a position that may or may not be CUTPOINT_UPPER or CUTPOINT_LOWER, determine whether this
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
update_cutpoint_virtual_atoms_if_connected( core::conformation::Conformation & conformation, core::Size const cutpoint_res, bool recurse = true );

void
show_atom_tree( kinematics::tree::Atom const & atom, Conformation const & conf, std::ostream & os );

void
replace_conformation_residue_copying_existing_coordinates(
	conformation::Conformation & conformation,
	Size const seqpos,
	chemical::ResidueType const & new_rsd_type
);

/// @brief Construct a variant of an existing conformation residue.
void add_variant_type_to_conformation_residue(
	conformation::Conformation & conformation,
	chemical::VariantType const variant_type,
	Size const seqpos );

/// @brief Construct a non-variant of an existing conformation residue.
void remove_variant_type_from_conformation_residue(
	conformation::Conformation & conformation,
	chemical::VariantType const variant_type,
	Size const seqpos );


void
add_lower_terminus_type_to_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
);


void
add_upper_terminus_type_to_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
);


void
remove_lower_terminus_type_from_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
);


void
remove_upper_terminus_type_from_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
);

void
build_tree(
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
);

/// @brief build a sub atom-tree for a jump edge and attach it to main atom-tree
void
build_jump_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
);

/// @brief build a sub atom-tree for a polymer edge and attach it to main atom-tree
void
build_polymer_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
);

/// @brief build a sub atom-tree for a chemical edge and attach it to main atom-tree
void
build_chemical_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
);


/// @brief  build the tree of atoms for this residue, anchored at root_atomno
void
build_residue_tree(
	int const root_atomno,
	conformation::Residue const & rsd,
	kinematics::AtomPointer1D & atom_ptr,
	bool const root_atom_is_jump_atom
	//bool const keep_1st_child_position = false
);

/// @brief  build_residue_tree function that uses the foldtree info
/// @brief  also used in build_tree to build the residue tree for the root residue
void
build_residue_tree(
	conformation::ResidueCOPs const & residues,
	conformation::Residue const & rsd,
	kinematics::FoldTree const & fold_tree,
	kinematics::AtomPointer1D & atom_ptr
);

/// @brief  Helper function for conformation routines
void
replace_residue_in_atom_tree(
	conformation::Residue const & new_rsd,
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
);

/// @brief  Inserts/ appends new residue subtree into an existing atomtree
/// @note  The foldtree must already have been changed to reflect the new residue
/// @note  The sequence position of the new residue is deduced from new_rsd.seqpos()
/// @note  This function handles renumbering of the atomtree if necessary
void
insert_residue_into_atom_tree(
	conformation::Residue const & new_rsd,
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
);

int
get_root_atomno(
	conformation::Residue const & rsd,
	int const dir // +1, -1, or "dir_jump"
);

Size
get_root_residue_root_atomno(
	conformation::Residue const & rsd,
	kinematics::FoldTree const & fold_tree
);

/// @brief  Get the atom-index of the atom to which the residue at position seqpos should be anchored.
int
get_anchor_atomno(
	conformation::Residue const & anchor_rsd,
	Size const seqpos,
	kinematics::FoldTree const & fold_tree
);


/// @brief get anchor atom to which the atom-tree of next residue in the edge is attached.
int
get_anchor_atomno(
	conformation::Residue const & rsd,
	int const dir // forward(1), backward(-1), or "dir_jump"
);

/// @brief  Use this routine to deduce atom indices of connect atoms in the tree
void
get_anchor_and_root_atoms(
	conformation::Residue const & anchor_rsd,
	conformation::Residue const & root_rsd,
	kinematics::Edge const & edge,
	Size & anchor_atomno,
	Size & root_atomno
);

/// @brief  Moves the first same-residue child of the jump atom corresponding to edge into first place in the child list
void
promote_sameresidue_child_of_jump_atom(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D const & atom_pointer
);


/// @brief  Moves the first same-residue child of the jump atom corresponding to edge into first place in the child list
void
promote_sameresidue_child_of_jump_atom(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
);

void
get_chemical_root_and_anchor_atomnos(
	conformation::Residue const & rsd_anchor,
	conformation::Residue const & rsd_root,
	Size & anchor_atom_no,
	Size & root_atom_no
);

/// @brief set up a map to match mainchain atoms from residue1 to residue2
void
setup_corresponding_atoms(
	id::AtomID_Map< id::AtomID > & atom_map,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
);

/// @brief Switch the disulfide state of a disulfide-forming residue (e.g. CYS->CYD or CYD->CYS or
/// DCYD->DCYS or DCYS->DCYD or whatnot).
/// @param[in] index Position of the residue to replace.
/// @param[in] cys_type_name3 The 3-letter name of the cys type to use: either CYS
///  or CYD.  DEPRECATED and kept only for backward-compatibility.
/// @param[inout] conf The conformation to modify
/// @return true if the replacement was successful, false otherwise.
bool change_cys_state(
	Size const index,
	std::string const &cys_type_name3,   //DEPRECATED and kept only for backward-compatibility.
	core::conformation::Conformation & conf
);

/// @brief Find whether there is a disulfide defined between two residues
bool is_disulfide_bond( core::conformation::Conformation const& conformation, core::Size residueA_pos, core::Size residueB_pos);

/// @brief Generate a list of all disulfide bonds in the conformation
void disulfide_bonds( core::conformation::Conformation const& conformation, utility::vector1< std::pair<core::Size,core::Size> > & disulfides );

/// @brief Gets a disulfide-forming residue's partner in the disulfide bond.
/// @author Vikram K. Mulligan, Baker lab (vmullig@uw.edu).
core::Size get_disulf_partner (
	core::conformation::Conformation const &conformation,
	core::Size const res_index
);

/// @brief Breaks a disulfide bond.
/// @author Vikram K. Mulligan, Baker lab (vmullig@uw.edu).
void break_disulfide(
	core::conformation::Conformation &conformation,
	core::Size const res1,
	core::Size const res2
);

/// @brief Introduce cysteines at the specified location and define a disulfide bond between them.
/// @details Does not do the repacking & minimization required to place the disulfide correctly.
void form_disulfide(
	core::conformation::Conformation & conformation,
	core::Size lower_res,
	core::Size upper_res,
	bool const preserve_d_residues=false,
	bool const force_d_residues=false
);

/// @brief Helper function for the form_disulfide function.
/// @details This function ensures that as a residue is mutated to a disulfide-bonding residue type,
/// all other variant types are preserved; it is used to avoid code duplication.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
void form_disulfide_helper(
	core::conformation::Conformation &conformation,
	core::Size const lower_res,
	core::chemical::ResidueTypeSetCOP restype_set,
	bool const preserve_d_residues,
	bool const force_d_residues
);

/// @brief Another helper function for the form_disulfide function.
/// @details Returns true if and only if the conformation is symmetric and upper_res is a symmetric copy of lower_res.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
bool
upper_is_symm_equivalent_of_lower(
	core::conformation::Conformation const &conformation,
	core::Size const lower_res,
	core::Size const upper_res
);

id::NamedAtomID
atom_id_to_named_atom_id(
	id::AtomID const & atom_id,
	conformation::Residue const & rsd
);


id::AtomID
named_atom_id_to_atom_id(
	id::NamedAtomID const & atom_id,
	conformation::Residue const & rsd
);

id::NamedStubID
stub_id_to_named_stub_id(
	id::StubID const & stub_id,
	conformation::Residue const & rsd
);


core::id::TorsionID find_bond_torsion_with_nearest_orientation(
	core::conformation::Conformation const & conf,
	utility::vector1< core::id::TorsionID > const & torsions,
	core::id::TorsionID const & query_torsion );


// Ring-related Functions /////////////////////////////////////////////////////

/// @brief  What is the attachment position of the query atom on the given ring?
core::uint position_of_atom_on_ring(
	Residue const & residue,
	core::uint query_atom,
	utility::vector1< core::uint > const & ring_atoms );

/// @brief  Is the query atom in this residue axial or equatorial to the given ring or neither?
chemical::rings::AxEqDesignation is_atom_axial_or_equatorial_to_ring(
	Residue const & residue,
	core::uint query_atom,
	utility::vector1< core::uint > const & ring_atoms );

/// @brief  Is the query atom in this residue axial or equatorial or neither?
//chemical::rings::AxEqDesignation is_atom_axial_or_equatorial( Residue const & residue, core::uint query_atom );

/// @brief Return the appropritate ResidueType for the virtual residue for the
/// "mode" (fullatom, centroid ...) the conformation is in.
///
/// When at all possible, use core::pose::virtual_type_for_pose() instead,
/// as that can use more pose-specific residue type information, if any.
chemical::ResidueTypeCOP
virtual_type_for_conf( core::conformation::Conformation const &conformation );

/// @brief Return the appropritate ResidueType for the inverse virtual residue for the
/// "mode" (fullatom, centroid ...) the conformation is in.
///
/// When at all possible, use core::pose::get_restype_for_pose() instead,
/// as that can use more pose-specific residue type information, if any.
chemical::ResidueTypeCOP
inv_virtual_type_for_conf( core::conformation::Conformation const &conformation );

/// @brief Get the center of the Residue
///
/// This computes an equally-weighted, all-atom (including virtuals and hydrogens) center
core::Vector
all_atom_center(
	core::conformation::Residue const & residue
);

/// @brief Returns a new residue based on a name
core::conformation::ResidueOP
get_residue_from_name(
	std::string const & name,
	std::string const & residue_type_set = "fa_standard" );

/// @brief Returns a new residue based on a name1
core::conformation::ResidueOP
get_residue_from_name1(
	char name1,
	bool is_lower_terminus = false,
	bool is_upper_terminus = false,
	bool d_aa = false,
	std::string const & residue_type_set = "fa_standard" );

/// @brief Given a residue and a connection id, get the heavyatom adjacent to the atom that makes that connection.
/// @details Chooses mainchain over non-mainchain, and heavyatoms over non-heavyatoms.  Returns true for FAILURE.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
get_second_atom_from_connection( core::Size & resno, core::Size & atomno, Residue const &rsd, Conformation const &conformation, core::Size const &conn_id );


/// @brief Get stub from orient_atoms
/// @details This defaults to a stub with CA, N, C for canonical amino acids
core::kinematics::Stub
get_stub_from_residue( core::conformation::Residue const & res );

/// @brief Get stub mapping that maps residue 1 to residue 2
/// @details Coordinate frames for each defined by orient_atoms which defaults to
//            CA, N, C for canonical amino acids
core::kinematics::RT
get_rt_from_residue_pair(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
);

} // conformation
} // core

#endif  // INCLUDED_core_conformation_util_HH
