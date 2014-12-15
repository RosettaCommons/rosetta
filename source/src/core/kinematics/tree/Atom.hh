// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/tree/Atom.hh
/// @brief  Kinematics Atom interface class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_tree_Atom_hh
#define INCLUDED_core_kinematics_tree_Atom_hh


// Unit headers
#include <core/kinematics/tree/Atom.fwd.hh>

// Package headers
// AUTO-REMOVED #include <core/kinematics/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>

// Numeric headers
#include <numeric/xyzMatrix.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector0.hh>
// AUTO-REMOVED #include <utility/vector1.hh> // DOH! switch all to vector1
#include <utility/pointer/ReferenceCount.hh>

#include <core/id/types.hh>
#include <utility/vector0_bool.hh>


namespace core {
namespace kinematics {
namespace tree {

/// Kinematics Atom interface class
class Atom : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Atom >
{

public: // Types


	typedef  PointPosition  Position;
	typedef  utility::vector0< AtomOP >  Atoms;
	typedef  Atoms::ConstIterator  Atoms_ConstIterator;
	typedef  Atoms::Iterator  Atoms_Iterator;

	typedef  numeric::xyzMatrix< Real > Matrix;


	// ids
	typedef  id::DOF_Type DOF_Type;
	typedef  id::DOF_ID DOF_ID;
	typedef  id::AtomID AtomID;
	typedef  id::AtomID_Mask AtomID_Mask;
	typedef  id::DOF_ID_Mask DOF_ID_Mask;


	// Types to prevent compile failure when std::distance is in scope
	typedef  void  iterator_category;
	typedef  void  difference_type;

protected: // Creation


	/// @brief Default constructor
	inline
	Atom()
	{}


	/// @brief Copy constructor
	inline
	Atom( Atom const & /*atom*/ ) : // PBHACK!!!
		ReferenceCount(),
		utility::pointer::enable_shared_from_this< Atom >()
	{}

public: // Creation


	/// @brief Destructor
	virtual
	~Atom()
	{}


public: // self pointers

	inline AtomCOP get_self_ptr() const { return shared_from_this(); }
	inline AtomOP get_self_ptr() { return shared_from_this(); }
	inline AtomCAP get_self_weak_ptr() const { return AtomCAP( shared_from_this() ); }
	inline AtomAP get_self_weak_ptr() { return AtomAP( shared_from_this() ); }


protected: // Assignment


	/// @brief Copy assignment
	inline
	Atom &
	operator =( Atom const & )
	{
		return *this;
	}


public: // Methods

	/// @brief Perform a depth-first traversal of the tree that would be effected by
	/// a DOF change from this atom.  Stop at atoms that have already been traversed.
	virtual
	void
	dfs(
		AtomDOFChangeSet & changeset,
		ResidueCoordinateChangeList & res_change_list,
		Size const start_atom_index
	) const = 0;


	/// @brief The atom must retrieve an appropriate stub from its parent; it is the root
	/// of the subtree being refolded
	virtual
	void
	update_xyz_coords() = 0;

	///////////////////////////////////////////////////////////////////////////
	// go back and forth between internal coords (DOF's)  and xyz coords
	/// @brief update xyz coords from stub and internal coords and
	virtual
	void
	update_xyz_coords(
		Stub & stub
	) = 0;

	/// @brief update internal coords from stub and xyz coords
	virtual
	void
	update_internal_coords(
		Stub & stub,
		bool const recursive = true
	) = 0;


	/// @brief calculate my input_stub from the current xyz's and use that input_stub to update my torsions
	virtual
	void
	update_internal_coords(
		bool const recursive
	) = 0;


	/// @brief update the stub without actually updating coordinates
	//Useful helper function for manipulating stubs
	virtual
	void
	update_stub(
		Stub & stub
	) const = 0;


	///////////////////////////////////////////////////////////////////////////
	/// @brief copy DOFs and xyz coords from src Atom
	virtual
	void
	copy_coords( Atom const & src ) = 0;


	///////////////////////////////////////////////////////////////////////////
	// access DOFs


	/// @brief get dof
	virtual
	Real
	dof(
		DOF_Type const type
	) const = 0;


	/// @brief set dof, use "set_" syntax since we have multiple dof's
	virtual
	void
	set_dof(
		DOF_Type const type,
		Real const value
	) = 0;


	/// @brief set dof, use "set_" syntax since we have multiple dof's -- for use in output-sensitive refold routine
	virtual
	void
	set_dof(
		DOF_Type const type,
		Real const value,
		AtomDOFChangeSet & set
	) = 0;


	/// @brief get Jump
	virtual
	Jump const &
	jump() const = 0;


	/// @brief set Jump
	virtual
	void
	jump(
		Jump const & jump_in
	) = 0;


	/// @brief set Jump -- for use in output-sensitive refolding
	virtual
	void
	jump(
		Jump const & jump_in,
		AtomDOFChangeSet & set
	) = 0;


	/// @brief copy atom with new memory allocation
	virtual
	AtomOP
	clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const = 0;


	///////////////////////////////////////////////////////////////////////////
	// for minimizing
	virtual
	void
	setup_min_map(
		DOF_ID & last_torsion,
		DOF_ID_Mask const & move_map,
		MinimizerMapBase & min_map
	) const = 0;


	virtual
	void
	get_dof_axis_and_end_pos(
		Vector & axis,
		Position & end_pos,
		DOF_Type const type
	) const = 0;


	///////////////////////////////////////////////////////////////////////////
	// miscellaneous inspection
	/// @brief atom is a jump atom?
	virtual
	bool
	is_jump() const = 0;

	/// @brief when other atoms are inserted insert after 1st child if available.
	/// --> this enables us to keep a stub of Downstream Jump atoms inside a single residue
	virtual
	bool
	keep_1st_child_pos() const = 0;


	/// @brief DoF should be fixed for this atom?
	///
	/// @details for DoFs that must be kept fixed due to topology of tree
	/// e.g., phi of stub_atoms for jump_atoms
	inline
	virtual
	bool
	keep_dof_fixed(
		DOF_Type const  //type
	) const
	{
		return false;
	}

	/// @brief dump out AtomID for this atom, its parent and all its offspring
	virtual
	void
	show() const = 0;

	/// @brief dump out AtomID for this atom, its parent and all its offspring up to n_level
	virtual
	void
	show(int const &) const = 0;

	/// @brief dihedral angle between two bonded children to this atom
	virtual
	Real
	dihedral_between_bonded_children(
		Atom const & child1,
		Atom const & child2
	) const = 0;


	///////////////////////////////////////////////////////////////////////////
	// update domain map
	virtual
	void
	update_domain_map(
		int & current_color,
		int & biggest_color,
		DomainMap & domain_map,
		AtomID_Mask const & dof_moved,
		AtomID_Mask const & atom_moved
	) const = 0;


	///////////////////////////////////////////////////////////////////////////
	// manage atom_list

	virtual
	Atoms_ConstIterator
	atoms_begin() const = 0;


	virtual
	Atoms_ConstIterator
	atoms_end() const = 0;


	virtual
	Atoms_Iterator
	atoms_begin() = 0;


	virtual
	Atoms_Iterator
	atoms_end() = 0;


	virtual
	Size
	n_atom() const = 0;


	// adds to the end, modulo the rule (which applies to the other methods
	// as well) that all the JumpAtoms come before the BondedAtoms
	virtual
	void
	append_atom( AtomOP ) = 0;


	virtual
	void
	delete_atom( AtomOP ) = 0;


	// inserts at the beginning
	virtual
	void
	insert_atom( AtomOP ) = 0;


	// tries to insert at the position specified by the second argument
	virtual
	void
	insert_atom( AtomOP, int const /*index*/ ) = 0;


	virtual
	void
	replace_atom(
		AtomOP const old_atom,
		AtomOP const new_atom
	) = 0;


	virtual
	AtomCOP
	get_nonjump_atom(
		Size const i
	) const = 0;


	virtual
	Size
	n_children() const = 0;


	virtual
	Size
	n_nonjump_children() const = 0;


	virtual
	AtomCOP
	child( Size const k ) const = 0;


	virtual
	AtomOP
	child( Size const k ) = 0;


	/// @brief the atom-index of this child
	virtual
	Size
	child_index( AtomCOP child ) const = 0;

	/// @brief the atom-index of this child
	virtual
	Size
	raw_child_index( Atom const * child ) const = 0;

	virtual
	bool
	downstream( AtomCOP atom1 ) const = 0;


public: // Properties


	/// @brief Atom identifier
	virtual
	AtomID const &
	id() const = 0;


	/// @brief AtomID assignment
	virtual
	void
	id( AtomID const & id_in ) = 0;


	/// @brief Atom identifier
	virtual
	AtomID const &
	atom_id() const = 0;


	/// @brief Position
	virtual
	Position const &
	position() const = 0;


	/// @brief Position assignment
	virtual
	void
	position( Position const & position_a ) = 0;


	/// @brief Position
	virtual
	Position const &
	xyz() const = 0;


	/// @brief Position assignment
	virtual
	void
	xyz( Position const & position_a ) = 0;


	/// @brief x coordinate
	virtual
	Length const &
	x() const = 0;


	/// @brief y coordinate
	virtual
	Length const &
	y() const = 0;


	/// @brief z coordinate
	virtual
	Length const &
	z() const = 0;


	/// @brief Distance to an Atom
	virtual
	Length
	distance( Atom const & atom ) const = 0;


	/// @brief Distance squared to an Atom
	virtual
	Length
	distance_squared( Atom const & atom ) const = 0;


	/// @brief Distance between two Atoms
	friend
	inline
	Length
	distance( Atom const & atom1, Atom const & atom2 );


	/// @brief Distance squared between two Atoms
	friend
	inline
	Length
	distance_squared( Atom const & atom1, Atom const & atom2 );


	/// @brief Transform atom and children by linear transformation
	virtual
	void
	transform_Ax_plus_b_recursive( Matrix const & A, Vector const & b, ResidueCoordinateChangeList & res_change_list ) = 0;

	/// @brief Parent atom pointer, NULL for root atom
	virtual
	AtomCOP
	parent() const = 0;

	virtual
	void
	get_path_from_root( utility::vector1< AtomCAP > & path ) const = 0;

	virtual
	bool
	atom_is_on_path_from_root( AtomCOP atm ) const = 0;

	/// @brief parent assignment
	virtual
	void
	parent( AtomAP parent_in ) = 0;

	/// @brief Parent atom pointer, NULL for root atom
	virtual
	AtomOP
	parent() = 0;

	/// @brief Get stub information
	virtual
	Stub
	get_stub() const = 0;

	virtual
	Stub
	get_input_stub() const = 0;

	virtual
	AtomCOP
	stub_atom1() const = 0;

	virtual
	AtomCOP
	stub_atom2() const = 0;

	virtual
	AtomCOP
	stub_atom3() const = 0;

	virtual
	AtomID const &
	stub_atom1_id() const = 0;

	virtual
	AtomID const &
	stub_atom2_id() const = 0;

	virtual
	AtomID const &
	stub_atom3_id() const = 0;

	virtual
	AtomCOP
	input_stub_atom0() const = 0;

	virtual
	AtomCOP
	input_stub_atom1() const = 0;

	virtual
	AtomCOP
	input_stub_atom2() const = 0;

	virtual
	AtomCOP
	input_stub_atom3() const = 0;

	virtual
	AtomID const &
	input_stub_atom0_id() const = 0;

	virtual
	AtomID const &
	input_stub_atom1_id() const = 0;

	virtual
	AtomID const &
	input_stub_atom2_id() const = 0;

	virtual
	AtomID const &
	input_stub_atom3_id() const = 0;

	// routines for navigating the tree
	virtual
	AtomCOP
	previous_sibling() const = 0;

	virtual
	AtomCOP
	previous_child(
		AtomCOP child
	) const = 0;

	virtual
	AtomOP
	next_child(
		AtomCOP child
	) = 0;

public:

	// The "raw" functions below are meant to be used only by AtomTree and only in performance critical
	// code where one is merely looking up data from the AtomTree, where the cost to increment the
	// reference count for each smart pointer in its constructor becomes prohibitive.  For example, the
	// AtomTree::torsion_angle function invokes stub_atom1, stub_atom2, and stub_atom3 creating dozens
	// of smart pointers.  This function is used following refold to copy data from the AtomTree into
	// the Residues; it suddenly started to take up a large percentage of the minimizer following the
	// move to smart pointers.  Why?  Atomic operations are much more expensive (which are what underlie
	// the boost::shared_ptr reference count data member) than simply incrementing and decrementing an
	// integer.  As a result, code that had not shown up in the profiler before now occupies much more
	// time than it used to.  It would be a misuse of these following functions to request and then hold
	// on to a raw pointer they returned.  They are meant for the lifetime of short, constant operations in
	// class AtomTree where nothing about the trees topology is changing.

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the parent atom pointer
	virtual
	Atom const *
	raw_parent() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the previous sibling pointer,
	/// i.e. the first child in the parent's children list to precede this atom.
	virtual
	Atom const *
	raw_previous_sibling() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the previous child pointer;
	virtual
	Atom const *
	raw_previous_child(
		Atom const * child
	) const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the fist stub atom
	virtual
	Atom const *
	raw_stub_atom1() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the second stub atom
	virtual
	Atom const *
	raw_stub_atom2() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the third stub atom
	virtual
	Atom const *
	raw_stub_atom3() const = 0;


	/// @brief Rapid (increment-of-reference-count-avoiding) access to the 0th input stub atom;
	virtual
	Atom const *
	raw_input_stub_atom0() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the 1st input stub atom;
	virtual
	Atom const *
	raw_input_stub_atom1() const = 0;

	/// @brief Rapid (increment-of-reference-count-avoiding) access to the 2nd input stub atom;
	virtual
	Atom const *
	raw_input_stub_atom2() const = 0;


	/// @brief Rapid (increment-of-reference-count-avoiding) access to the 3rd input stub atom;
	virtual
	Atom const *
	raw_input_stub_atom3() const = 0;


	/// @brief Rapid (increment-of-reference-count-avoiding) access to the ith non-jump atom in this atom's
	/// list of children.
	virtual
	Atom const *
	raw_get_nonjump_atom(
		Size const i
	) const = 0;


	virtual
	bool
	stub_defined() const = 0;


protected: // Methods


	// when subtrees have changed their coordinates
	virtual
	void
	update_child_torsions(
		AtomOP const child
	) = 0;


	virtual
	Atoms_ConstIterator
	nonjump_atoms_begin() const = 0;


	virtual
	Atoms_Iterator
	nonjump_atoms_begin() = 0;

}; // Atom


/// @brief Distance between two Atoms
inline
Length
distance( Atom const & atom1, Atom const & atom2 )
{
	return atom1.distance( atom2 );
}

/// @brief Distance squared between two Atoms
inline
Length
distance_squared( Atom const & atom1, Atom const & atom2 )
{
	return atom1.distance_squared( atom2 );
}



} // namespace tree
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_Atom_HH
