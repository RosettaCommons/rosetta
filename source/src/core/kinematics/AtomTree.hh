// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/AtomTree.hh
/// @brief  Atom tree class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_AtomTree_hh
#define INCLUDED_core_kinematics_AtomTree_hh


// Unit headers
#include <core/kinematics/AtomTree.fwd.hh>

// Package headers
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.hh> // apl temp, until all of AtomTree's methods are moved into its .cc file
#ifdef WIN32
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomWithDOFChange.hh>
#endif

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>

// Numeric headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <map>


namespace core {
namespace kinematics {


/// @brief The AtomTree class is a wrapper for a tree of kinematic Atoms.
///
/// @details The responsibilities of the class are:
///
/// @li 1. Maintain a map from AtomIDs to Atoms for fast lookup.
///
/// @li 2. Keep the internal and xyz coords of the Atoms in sync.
///    Note that this sync'ing is handled
///    in a lazy fashion, analogous to the way the current pose handles
///    refolding. As a result, getting and setting of coords can trigger
///    coordinate updates.
///
/// @li 3. Gatekeep modification of internal or xyz coords for the Atoms
///    (necessary for #2).
///
/// @li 4. Serve out the xyz coords for updating a Pose/Conformation object.
///
/// @li 5. Identify DOF_IDs that correspond to torsion angles (e.g., phi, psi, chi, nu)
///    specified by 4 AtomIDs. Do this in a way that enables fast
///    lookup, e.g., caching the results of previous calls (?).
///    Note that some torsions will not correspond exactly to a DOF of
///    an Atom; e.g., chi1 when we are folding c->n, necessitating
///    calculation of an offset. We should be able to handle getting/setting
///    of these torsion angles as well as handing back the DOF_ID
///    (the last one is necessary when setting up the DOF_IDMask for
///    minimization given a MoveMap object).
///
/// See @ref atomtree_overview "AtomTree overview and concepts" for details.
class AtomTree : public utility::pointer::ReferenceCount
{
public: // Types
	// ids
	typedef  id::AtomID AtomID;
	typedef  id::AtomID_Mask AtomID_Mask;
	typedef  id::DOF_Type DOF_Type;
	typedef  id::DOF_ID DOF_ID;
	typedef  id::DOF_ID_Mask DOF_ID_Mask;
	typedef  id::StubID StubID;

	typedef std::vector< AtomID > AtomIDs;

	// for fragment insertions
	typedef std::map< AtomID, Vector > FragXYZ;
	typedef std::map< StubID, RT > FragRT;

	typedef tree::Atom    Atom;
	typedef tree::AtomOP  AtomOP;
	typedef tree::AtomCOP AtomCOP;

public: // Creation

	////////////////////////////
	// note that we steal the atoms, i.e., they're incorporated into the AtomTree,
	// not cloned
	/// @brief construction: from a tree of atoms;
	AtomTree( AtomPointer2D const & new_atom_pointer, bool const from_xyz = true );

	/// @brief default constructor
	AtomTree();

	/// @brief Destructor
	virtual
	~AtomTree();

	/// @brief copy constructor
	AtomTree( AtomTree const & src );

public: // Assignment

	/// @brief Copy assignment, makes complete copy of another AtomTree
	AtomTree &
	operator=( AtomTree const & src );

public:
	/// @brief Weak-pointer setter.  The object that instantiates an owning pointer to an AtomTree object
	/// must hand that AtomTree a weak pointer to itself so that the AtomTree may share that weak pointer
	/// with other AtomTrees.  Such sharing allows for crucial speedups when copying between AtomTrees.
	/// If the object that instantiates this AtomTree does not provide it with a pointer-to-self, the
	/// AtomTree will still function, but it will not share its pointers properly.
	void set_weak_pointer_to_self( AtomTreeCAP self_pointer );

public: // Methods

	/// @brief number of residues
	Size
	size() const
	{
		return atom_pointer_.size();
	}

	/// the atom with this AtomID in the current AtomID_Map?
	/// @brief true only if AtomID is in the range of the map and its Atom pointer has value
	bool
	has( AtomID const & id ) const
	{
		return ( atom_pointer_.has( id ) && ( atom_pointer_[ id ] != 0 ) );
	}


	void
	replace_residue_subtree(
		id::BondID const & incoming,
		utility::vector1< id::BondID > const & outgoing,
		AtomPointer1D const & new_atoms // this will be the new slice of our atom_pointer
	);


	/// @brief  Useful for guaranteeing that a stub remains within a single residue
	void
	promote_sameresidue_nonjump_child( AtomID const & parent_atom_id );

	/// @brief  Deletes atoms for seqpos.
	/// Does not renumber other atoms -- need to call update_sequence_numbering for that
	/// designed for the simple case of 1 incoming connecting and 0 or 1 outgoing connection,
	/// where the desired behavior is to rewire the outgoing connection in place of seqpos' tree
	void
	delete_seqpos( Size const seqpos );

	/// @brief updates the Atom's AtomID's and the atom_pointer array
	void
	update_sequence_numbering(
		Size const new_size,
		utility::vector1< int > const & old2new
	);

	/// @brief replaces the entire tree
	void
	replace_tree( AtomPointer2D const & new_atom_pointer, bool const from_xyz = true );


	/// @brief copy the internal and xyz coords from src tree, which should have the same topology as me
	void
	copy_coords(
		AtomTree const & src
	);


	/// @brief set a specific DOF in the tree
	void
	set_dof( DOF_ID const & id, Real const setting );

	/// @brief set a specific atom xyz position
	void
	set_xyz( AtomID const & id, PointPosition const & xyz );

	/// @brief simultaniously set several atom xyz positions
	void
	batch_set_xyz( utility::vector1<AtomID> const & id, utility::vector1<PointPosition> const & xyz );

	/// @brief set a specific jump transform
	void
	set_jump( AtomID const & id, Jump const & jump );

	/// @brief set a specific jump transform and immediately refold downstream atoms
	void
	set_jump_now( AtomID const & id, Jump const & jump );

  /// @brief find the atom giving rise to the jump connecting two stubs.
  /// @param direction is set to 1 or -1 depending on the direction of the jump
  AtomID
  get_jump_atom_id(
    StubID const& stub_id1,
    StubID const& stub_id2,
    int& direction
  ) const;


	/// @brief  Set the transform between two stubs, returns the atomid of the jump atom which moved (for book-keeping)
	AtomID
	set_stub_transform(
		StubID const & stub_id1,
		StubID const & stub_id2,
		RT const & target_rt
	);

	/// @brief  get the transform between two stubs
	RT
	get_stub_transform(
		StubID const & stub_id1,
		StubID const & stub_id2
	) const;

	/// @brief set a torsion angle "setting" to a specifc DOF derived by the four atoms
	DOF_ID
	set_torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		Real const setting,
		bool const quiet=false
	);

	DOF_ID
	set_bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		Real const setting
	);

	DOF_ID
	set_bond_length(
		AtomID const & atom1,
		AtomID const & atom2,
		Real const setting
	);


	/// @brief generates a "domain_map" defining the rigid body regions
	/// whose internal coords have not changed, according to the
	/// informaiton in the two bool Mask's
	/// update domain map from dof_moved and xyz_moved
	void
	update_domain_map(
		DomainMap & domain_map,
		AtomID_Mask const & dof_moved,
		AtomID_Mask const & xyz_moved
	) const;


	/// @brief clear the content of an AtomTree object, delete root atom and all pointers
	void
	clear();

	void
	insert_fragment(
		StubID const & instub_id,
		FragRT const & outstub_transforms,
		FragXYZ const & frag_xyz,
		utility::vector1< AtomID > & moving_atoms
	);


	/// @brief The AtomTree provides to the Conformation object a list of residues
	/// whose xyz coordinates have changed.  When the Conformation has finished reading off
	/// residues that have changed from the AtomTree, and has copied the coordinates of
	/// those residues into its conformation::Residue objects, it informs the AtomTree
	/// to reset this list by a call to mark_changed_residues_registered
	ResidueListIterator
	residue_xyz_change_list_begin() const;

	ResidueListIterator
	residue_xyz_change_list_end() const;

	/// @brief The AtomTree provides a list of residues who's xyz coordinates have changed
	/// to the Conformation object.  When the Conformation has finished reading off residues
	/// that have changed from the AtomTree, and has copied the coordinates of those residues
	/// into its conformation::Residue objects, it informs the AtomTree to reset this list
	/// by a call to mark_changed_residues_registered
	void
	note_coordinate_change_registered() const;

public: // Properties

	/// @brief is there any atom in the tree yet?
	bool
	empty() const
	{
		return ( root_ == 0 );
	}

	/// @brief  const-access to the root of the tree
	AtomCOP
	root() const
	{
		return root_;
	}


	// accessors -- which may trigger coordinate updates
	/// @brief get value of DOF( PHI, THETA, D, RB1, ....)  given its DOF_ID
	Real
	dof( DOF_ID const & id ) const;

	/// @brief get xyz position of an atom given its AtomID
	PointPosition const &
	xyz( AtomID const & id ) const;


	Atom const &
	atom( AtomID const & id ) const;

	Atom const &
	atom_dont_do_update( AtomID const & id ) const;


	Jump const &
	jump( AtomID const & id ) const;

	/// @brief a wrapper function to get the DOF_ID of a torsion angle given those four atoms which define this torsion
	/// @details Another version of this function also calculates an offset value, which is not needed here.
	DOF_ID
	torsion_angle_dof_id(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		bool const quiet=false
	) const
	{
		Real offset;
		return torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset, quiet );
	}

	/// @brief get the DOF_ID of a torsion angle given those four atoms which define this torsion
	DOF_ID
	torsion_angle_dof_id(
		AtomID const & atom1_in_id,
		AtomID const & atom2_in_id,
		AtomID const & atom3_in_id,
		AtomID const & atom4_in_id,
		Real & offset,
		bool const quiet=false
	) const;


	/// @brief get the DOF_ID of a bond angle given those 3 atoms which define this torsion
	DOF_ID
	bond_angle_dof_id(
		AtomID const & atom1_in_id,
		AtomID const & atom2_in_id,
		AtomID const & atom3_in_id,
		Real & offset
	) const;


	/// @brief get the DOF_ID of a bond length given those 2 atoms which define this torsion
	DOF_ID
	bond_length_dof_id(
		AtomID const & atom1_in_id,
		AtomID const & atom2_in_id
	) const;

	/// @brief calculate torsion angle defined by four atoms in the atom tree
	Real
	torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4
	) const;

	Real
	bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3
	) const;

	Real
	bond_length(
		AtomID const & atom1,
		AtomID const & atom2
	) const;


	void
	set_jump_atom_stub_id( StubID const & id );

	Stub
	stub_from_id( StubID const & id ) const
	{
		if ( id.center().valid() ) {
			return Stub( xyz( id.center() ), xyz( id.atom1 ), xyz( id.atom2 ), xyz( id.atom3 ) );
		} else {
			return Stub( xyz( id.atom1 ), xyz( id.atom2 ), xyz( id.atom3 ) );
		}
	}

private:

	/// @brief When an atom tree copies the topology of another atom tree, it must
	/// register itself as a topological observer of that other tree.  When the other
	/// tree changes its topology, then the tree being observed must notify its
	/// observers that they are no longer a topological copy of this tree.  An atom
	/// tree may only be the topological copy of a single other atom tree, though several
	/// atom trees may be copies of a single atom tree.
	void attach_topological_observer( AtomTreeCAP observer ) const;

	/// @brief When an atom tree changes its topology, it must inform all of its
	/// observers that they are no longer the same topology as this tree.
	void notify_topological_change( AtomTreeCAP observee ) const;

	/// @brief When an atom tree observing this tree decides it wants to become an observer
	/// of another tree, it must notify the tree that it formerly observed of this change.
	void detatch_topological_observer( AtomTreeCAP observer ) const;

public:
	/// Functions only necessary for unit tests

	/// @brief For testing purposes only: report the address of the AtomTree this tree
	/// is a topological copy of.  The fact that AtomTrees keep track of other atom trees
	/// is "private" in the sense that no other class needs to worry about it.  However,
	/// to *test* that the topological match algorithm is working properly, this private
	/// data needs to be readable.  Do not use this function outside of the unit tests.
	AtomTreeCAP
	topological_match_to() const {
		return topological_match_to_;
	}

	/// @brief For testing purposes only: report the list of observer AtomTrees that
	/// are topological copies of this tree.  The fact that AtomTrees keep track of
	/// other atom trees is "private" in the sense that no other class needs to worry
	/// about it.  However, to *test* that the topological match algorithm is working
	/// properly, this private data needs to be readable.  Do not use this function
	/// outside of the unit tests.
	utility::vector1< AtomTreeCAP > const &
	topological_observers() const {
		return topological_observers_;
	}


private: // Helper Methods for fragment insertion


	/// @brief  Deduce root_ from atom_pointer_ -- look for atom with atom->parent() == 0
	void
	find_root_from_atom_pointer();

	void
	get_frag_atoms(
		StubID const & id,
		FragXYZ const & frag_xyz,
		AtomCOP & frag_atom,
		AtomCOP & nonfrag_atom
	) const;


	StubID
	get_frag_pseudo_stub_id(
		AtomID const & id,
		FragXYZ const & frag_xyz,
		bool & fail
	) const;


	Stub
	get_frag_local_stub(
		StubID const & stubid,
		FragXYZ const & frag_xyz,
		bool & fail
	) const;

	Vector
	get_frag_local_xyz(
		AtomID const & id,
		FragXYZ const & frag_xyz,
		bool & fail
	) const;


	Vector
	get_frag_descendant_local_xyz(
		AtomCOP atom,
		FragXYZ const & frag_xyz,
		bool & fail
	) const;


	Vector
	get_frag_parent_local_xyz(
		AtomCOP child,
		FragXYZ const & frag_xyz,
		bool & fail
	) const;


	void
	insert_single_fragment(
		StubID const & instub_id,
		FragRT const & outstub_transforms,
		FragXYZ const & frag_xyz,
		utility::vector1< AtomID > & moving_atoms
	);

private: // Methods

	/// @brief bookkeeping -- set the Atoms' atomIDs from the atom_pointer_ array
	void
	update_atom_ids_from_atom_pointer();


	/// @brief  Convenience when we want an Atom*
	AtomOP
	atom_pointer( AtomID const & id )
	{
		return atom_pointer_[ id ]();
	}

	/// @brief  Convenience when we want an Atom*
	AtomCOP
	atom_pointer( AtomID const & id ) const
	{
		return atom_pointer_[ id ]();
	}


	// these two private functions are for maintaining synchrony between the internal and xyz coords

	/// @brief  Update the internal coordinates using the xyz (cartesian) coords
	void
	update_internal_coords() const;

	/// @brief  Update the xyz coordinates using the internal coords
	void
	update_xyz_coords() const;


	/// @brief  Notify self of new tree topology
	/// Useful if we move to caching some things that depend on the tree
	/// @note Probably need to go through and put more calls of this guy
	void
	set_new_topology();

private: // Fields

	/// @brief A weak pointer to self (this).
	/// @details The weak pointer must be provided to the AtomTree immediately after creation.
	AtomTreeCAP this_weak_ptr_;

	/// @brief Root Atom
	AtomOP root_;

	/// @brief Atom pointers map (map[AtomID] = AtomPointer)
	AtomPointer2D atom_pointer_;

	/// @brief List of the jump atom ID's, excluding the root. Order matters (for movemap indexing)
	//utility::vector1< AtomID > jump_atoms_; // NOT HERE YET

	/// @brief Internal coords out of date?
	mutable bool internal_coords_need_updating_;

	/// @brief XYZ coords out of date?
	mutable bool xyz_coords_need_updating_;

	/// @brief pointer to the atom tree this tree has an exact topological match to
	/// since that tree was the last tree copied from without subsequence topological
	/// modifications -- or at most one modification when that tree copied this
	/// tree's topology
	mutable AtomTreeCAP topological_match_to_;

	/// @brief pointers to all atom trees that are observing this tree's topology.
	/// On topological changes (including the destruction of this tree!),
	/// each of these trees have their topological_match_to_
	/// pointers set to null and this list is cleared.
	mutable utility::vector1< AtomTreeCAP > topological_observers_;

	/// @brief A list of the atoms that have had changed DOFs since the last refold.
	mutable AtomDOFChangeSet dof_changeset_;

	/// @brief A list of residues that have had xyz coordinate changes since the last
	/// time the owning Conformation object has asked for an update.
	ResidueCoordinateChangeListOP external_coordinate_residues_changed_;

	/// @ (SOON) A list of residues that have had DOF changes since the last
	/// time the owning Conformation object has asked for an update.
	//ResidueCoordinateChangeListOP internal_coordinate_residues_changed_;

}; // AtomTree

} // namespace kinematics
} // namespace core

#endif // INCLUDED_core_kinematics_AtomTree_HH
