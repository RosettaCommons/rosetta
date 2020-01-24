// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/tree/JumpAtom.hh
/// @brief  Jump atom
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_tree_JumpAtom_hh
#define INCLUDED_core_kinematics_tree_JumpAtom_hh


// Package headers
#include <core/kinematics/Jump.hh>
#include <core/kinematics/tree/Atom_.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {
namespace tree {

/// @brief an atom who are connected to its parent via rigid-body transformation "Jump"
///
/// See @ref atomtree_overview "AtomTree overview and concepts" for details.
///
class JumpAtom : public Atom_
{
private: // Types


	typedef  Atom_  Super;


public:

	/// @brief Perform a depth-first traversal of the tree that would be effected by
	/// a DOF change from this atom.  Stop at atoms that have already been traversed.
	void
	dfs(
		AtomDOFChangeSet & changeset,
		ResidueCoordinateChangeList & res_change_list,
		Size const start_atom_index
	) const override;

	///////////////////////////////////////////////////////////////////////////
	// go back and forth between DOFs and coords

	/// @brief The atom must retrieve an appropriate stub from its parent; it is the root
	/// of the subtree being refolded.  Valid only if this atom is the maximal root of a subtree
	/// requiring coordinate updates -- if any ancestor of this atom requires a coordinate update,
	/// then the Stub this atom generates for itself will be invalid.
	void
	update_xyz_coords() override;

	/// @brief update this atom's xyz position
	void
	update_xyz_coords(
		Stub & stub
	) override;

	using Atom_::update_internal_coords;

	/// update the jump info
	void
	update_internal_coords(
		Stub & stub,
		bool const recursive = true ) override;


	///////////////////////////////////////////////////////////////////////////
	// access DOFs

	/// @brief set a degree of freedom for jump
	void
	set_dof(
		DOF_Type const type,
		Real const value
	) override;

	/// @brief set degrees of freedom (internal coordinates).  For use in
	/// output-sensitive refold subroutine.
	void
	set_dof(
		DOF_Type const type,
		core::Real const value,
		AtomDOFChangeSet & changeset
	) override;


	/// @brief get a degree of freedom from jump
	Real
	dof(
		DOF_Type const type
	) const override;

	/// @brief access the jump
	Jump const &
	jump() const override;

	/// @brief set the jump
	void
	jump( Jump const & jump_in ) override;


	/// @brief set the jump.  For use with output-sensitive refold subroutine.
	void
	jump( Jump const & jump_in, AtomDOFChangeSet & changeset ) override;


	/// @brief copy this atom
	AtomOP
	clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const override;

	// get inversion of the wrapped jump
	std::pair<bool,bool>
	get_inversion() {
		return std::make_pair<bool,bool>(jump_.get_invert_upstream(), jump_.get_invert_downstream());
	}

	void
	steal_inversion(AtomOP steal_from) override;

	///////////////////////////////////////////////////////////////////////////
	/// @brief for minimizing, add DOF(RB) for a JumpAtom into the MinimizerMap
	void
	setup_min_map(
		DOF_ID & last_torsion,
		DOF_ID_Mask const & allow_move,
		MinimizerMapBase & min_map
	) const override;

	/// @brief get rotation axis and end_pos for a JumpAtom.
	void
	get_dof_axis_and_end_pos(
		Vector & axis,
		Position & end_pos,
		DOF_Type const type
	) const override;


	///////////////////////////////////////////////////////////////////////////
	// miscellaneous inspection

	/// @brief a jump atom is a jump? of course yes!!!
	inline
	bool
	is_jump() const override { return true; }

	///\brief when other atoms are inserted insert after 1st child if available.
	/// --> this enables us to keep a stub of Downstream Jump atoms inside a single residue
	bool
	keep_1st_child_pos() const override { return false; }

	/// whether a jump should be fixed in some special cases
	bool
	keep_dof_fixed(
		DOF_Type const type
	) const override;


	///////////////////////////////////////////////////////////////////////////

	/// @brief copy DOFs, xyz's
	void
	copy_coords( Atom const & src ) override;


	//protected:

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// protected methods
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	// useful helper function for manipulating stubs
	/// @brief  update the stub without actually updating coordinates
	/** since for a jump atom, update internal coords or xyz dont change input
	jump, so we do not do anything here*/
	inline
	void
	update_stub(
		Stub & //stub
	) const override
	{} // stub doesnt change


public: // Properties


	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom1 of a jump-atom
	/** it is itself if a stub can be defined for it. Otherwise it is parent*/
	inline
	AtomCOP
	stub_atom1() const override
	{
		return ( stub_defined() ? get_self_ptr() : parent() );
	}


	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom2 of a jump-atom
	/** it is its first bonded child if a stub can be defined for it. Otherwise
	it is parent's stub_atom2. */
	inline
	AtomCOP
	stub_atom2() const override
	{
		if ( stub_defined() ) {
			return get_nonjump_atom(0);
		}
		AtomCOP parent_op = parent();
		if ( parent_op ) {
			return parent_op->stub_atom2();
		}
		return nullptr;
	}
	/////////////////////////////////////////////////////////////////////////////
	/// @brief stub_atom3 of a jump atom
	/** it is its child's child or its second child if a stub can be defined for it,
	otherwise it is its parent's stub_atom3 */
	inline
	AtomCOP
	stub_atom3() const override
	{
		//std::cout << "stub_atom3: " << this << ' ' << parent_ << std::endl();
		if ( stub_defined() ) {
			AtomCOP first( get_nonjump_atom(0) );
			AtomCOP second( first->get_nonjump_atom(0) );
			if ( second ) {
				return second;
			} else {
				return get_nonjump_atom(1);
			}
		}
		AtomCOP parent_op = parent();
		if ( parent_op ) {
			return parent_op->stub_atom3();
		}
		return nullptr;
	}

public:

	Atom const *
	raw_stub_atom1() const override;

	Atom const *
	raw_stub_atom2() const override;

	Atom const *
	raw_stub_atom3() const override;

private: // Fields


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// data
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	/// @brief Jump
	/**
	A jump atom is connected to its parent via rigid-body transformation("jump").
	It requires two stubs to define a jump, one is the parent atoms's stub and
	the other is the stub centered at this jump atom, which requires at least
	three atoms on the jump atom's side (including itself). For example, a stub
	is defined from B-A-C or A-B-C (A is the jump atom, B and C are its offspring).
	If less than 3 atoms on the jump atom's side, i.e., stub_defined() == False,
	this atom will just use its parent stub.
	*/
	Jump jump_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // JumpAtom

typedef utility::pointer::shared_ptr< JumpAtom > JumpAtomOP;
typedef utility::pointer::shared_ptr< JumpAtom const > JumpAtomCOP;
typedef utility::pointer::weak_ptr< JumpAtom > JumpAtomAP;
typedef utility::pointer::weak_ptr< JumpAtom const > JumpAtomCAP;

} // namespace tree
} // namespace kinematics
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_kinematics_tree_JumpAtom )
#endif // SERIALIZATION


#endif // INCLUDED_core_kinematics_JumpAtom_HH
