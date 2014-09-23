// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/tree/BondedAtom.hh
/// @brief  Kinematics
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_tree_BondedAtom_hh
#define INCLUDED_core_kinematics_tree_BondedAtom_hh


// Package headers
#include <core/kinematics/tree/Atom_.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <cassert>


namespace core {
namespace kinematics {
namespace tree {

extern Jump BOGUS_JUMP; // for return statement to keep compiler happy

/// @brief an atom which are bonded to its parent, derived from Atom_
///
/// See @ref atomtree_overview "AtomTree overview and concepts" for details.
///
class BondedAtom : public Atom_
{

public:

	// default constructor
	BondedAtom() :
		phi_(0.0),
		theta_(0.0),
		d_(0.0),
		dof_change_propagates_to_younger_siblings_( false )
	{}

private: // Types


	typedef  Atom_  Super;


public:

	/// @brief Perform a depth-first traversal of the tree that would be effected by
	/// a DOF change from this atom.  Stop at atoms that have already been traversed.
	/// Will recurse on younger siblings if a phi on this atom has changed.
	virtual
	void
	dfs(
		AtomDOFChangeSet & changeset,
		ResidueCoordinateChangeList & res_change_list,
		Size const start_atom_index
	) const;

	///////////////////////////////////////////////////////////////////////////
	// go back and forth between DOFs and coords

	/// @brief The atom must retrieve an appropriate stub from its parent; it is the root
	/// of the subtree being refolded.  Valid only if this atom is the maximal root of a subtree
	/// requiring coordinate updates -- if any ancestor of this atom requires a coordinate update,
	/// then the Stub this atom generates for itself will be invalid.
	virtual
	void
	update_xyz_coords();


	/// @brief update cartesian coordinates for this atom from its input stub and internal cooridnates
	virtual
	void
	update_xyz_coords(
		Stub & stub
	);

	using Atom_::update_internal_coords;

	/// @brief update internal coordinates for this atom from its xyz position and input stub
	virtual
	void
	update_internal_coords(
		Stub & stub,
		bool const recursive = true
	);


	// useful helper function for manipulating stubs
	/// @brief update the stub without actually updating coordinates
	virtual
	void
	update_stub(
		Stub & stub
	) const;


	///////////////////////////////////////////////////////////////////////////
	// access DOFs

	/// @brief set degrees of freedom (internal coordinates)
	virtual
	void
	set_dof(
		DOF_Type const type,
		core::Real const value
	);

	/// @brief set degrees of freedom (internal coordinates). For use in
	/// output-sensitive refold subroutine.
	virtual
	void
	set_dof(
		DOF_Type const type,
		core::Real const value,
		AtomDOFChangeSet & changeset
	);

	/// @brief get degrees of freedom
	virtual
	core::Real
	dof(
		DOF_Type const type
	) const;

	/// @brief abort if attempt to get jump for a bonded atom
	inline
	virtual
	Jump const &
	jump() const { abort_bad_call(); return BOGUS_JUMP; /* we never get here */ }


	/// @brief abort if attempt to set jump for a bonded atom
	inline
	virtual
	void
	jump( Jump const & /* jump_in */ ) { abort_bad_call(); }

	/// @brief abort if attempt to set jump for a bonded atom
	inline
	virtual
	void
	jump( Jump const & /* jump_in */, AtomDOFChangeSet & /*changeset*/ ) { abort_bad_call(); }


	/// @brief copy this atom
	virtual
	AtomOP
	clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const;


	///////////////////////////////////////////////////////////////////////////
	///@brief  for minimizing,add DOF(PHI,THETA,D) for a BondedAtom into the MinimizerMap
	virtual
	void
	setup_min_map(
		DOF_ID & last_torsion,
		DOF_ID_Mask const & allow_move,
		MinimizerMapBase & min_map
	) const;

	///@brief get rotation axis and end_pos for a BondedAtom.
	virtual
	void
	get_dof_axis_and_end_pos(
		Vector & axis,
		Position & end_pos,
		DOF_Type const type
	) const;


	///////////////////////////////////////////////////////////////////////////
	// miscellaneous inspection

	/// @brief bonded atom is a jump? of course not!!!
	inline
	virtual
	bool
	is_jump() const { return false; }

	///\brief when other atoms are inserted insert after 1st child if available.
	/// --> this enables us to keep a stub of Downstream Jump atoms inside a single residue
	virtual
	bool
	keep_1st_child_pos() const { return false; }


	/// @brief whether a DOF for this atom should be fixed?
	virtual
	bool
	keep_dof_fixed(
		DOF_Type const type
	) const;


	///////////////////////////////////////////////////////////////////////////

	/// @brief copy DOFs, xyz's
	virtual
	void
	copy_coords( Atom const & src );


public: // Properties


	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom1 of a bonded atom
	/** it is itself */
	inline
	AtomCOP
	stub_atom1() const
	{
		return get_self_ptr();
	}


	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom2 of a bonded atom
	/** it is its parent */
	inline
	AtomCOP
	stub_atom2() const
	{
		return parent();
	}


	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom3 of a bonded atom
	/**
			- if this atom's parent is a bonded_atom it is this atom's parent's parent.
			- if this atom's parent is a jump atom, it is this atom's first non-jump
				sibling or its second non-jump sibling (if it itself is the first) or
				its first non-jump child (if it does not have any sibling)
	*/
	inline
	AtomCOP
	stub_atom3() const
	{
		//std::cout << "stub_atom3: " << this << ' ' << parent_ << std::endl();
		AtomCOP parent_op = parent(); // must have parent
		if ( parent_op->is_jump() ) {
			assert( parent_op->stub_defined() ); // weird behavior otherwise
			AtomCOP p_stub2( parent_op->stub_atom2() );
			AtomID const & p_stub2_id( p_stub2->id() );
			if ( id() == p_stub2_id ) {
				// very special case!!
				return parent_op->stub_atom3();
			} else {
				return p_stub2;
			}
		} else {
			return parent_op->stub_atom2();
		}
	}


private: // Fields


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// data
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	/// @brief DOF properties of a bonded atom
	/**
		 - a bonded atom is an atom who connects to its parent by a "bond" (covalent
		 or virtual) and therefore its position is defined by three internal
		 coordinates: d_, theta_, and phi_.
		 - d_ is the bond distance between this atom and its parent.
		 - theta_ is the bond angle between this atom(A), its parent's stub_atom1(B)
		 and its parent's stub_atom2(C), i.e., angle in between B->A ^ C->B.
		 - phi_ is either the torsion angle defined by A, B, C and D (C's parents),
		 or the improper angle defined by A, B, C and D (B's child and A's previous sibling).
	 */
	Real phi_, theta_, d_;

	/// @brief Track whether a dof change from this node (since the last update_xyz)
	/// induces a coordinate change for this node's younger siblings.
	bool dof_change_propagates_to_younger_siblings_;


}; // BondedAtom

typedef utility::pointer::owning_ptr< BondedAtom > BondedAtomOP;
typedef utility::pointer::owning_ptr< BondedAtom const > BondedAtomCOP;
typedef utility::pointer::access_ptr< BondedAtom > BondedAtomAP;
typedef utility::pointer::access_ptr< BondedAtom const > BondedAtomCAP;

} // namespace tree
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_BondedAtom_HH
