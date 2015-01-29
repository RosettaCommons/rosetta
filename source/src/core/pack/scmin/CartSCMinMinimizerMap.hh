// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/CartSCMinMinimizerMap.hh
/// @brief  Class for identifying the sidechain DOFs in the AtomTree which are free during
///         any particular call to the minimizer.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_CartSCMinMinimizerMap_hh
#define INCLUDED_core_pack_scmin_CartSCMinMinimizerMap_hh

// Unit headers
#include <core/pack/scmin/CartSCMinMinimizerMap.fwd.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>
#include <core/pack/scmin/CartSCMinMultifunc.hh>
#include <core/optimization/Multifunc.fwd.hh>

// Package Headers
#include <core/pack/scmin/AtomTreeCollection.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/optimization/types.hh>
#include <core/optimization/DOF_Node.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>


namespace core {
namespace pack {
namespace scmin {

/// @brief
class CartSCMinMinimizerMap : public SCMinMinimizerMap
{
public:
	typedef optimization::DOF_Node   DOF_Node;
	typedef optimization::DOF_NodeOP DOF_NodeOP;

public:
	CartSCMinMinimizerMap();
	virtual ~CartSCMinMinimizerMap();

	/// @brief the CartSCMinMinimizerMap has to know how many residues are in the pose; this allows
	/// it to do O(1) updates to its DomainMap -- this function costs O(N).
	virtual
	void set_total_residue( Size total_residue );

	/// @brief Disable the minimization for all residues.  Ammortized O(1).
	virtual
	void clear_active_dofs();

	/// @brief Activate all the dofs for a particular residue.  Ammortized O(1).
	virtual
	void activate_residue_dofs( Size resindex );

	/// @brief Invoked during the depth-first traversal through the AtomTree.  The AtomTree
	/// is indicating that a particular torsion is dependent on another torsion.  Record
	/// that fact.
	virtual
	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	);

	/// @brief Invoked during the depth-first traversal through the AtomTree; the atom
	/// tree is indicating that a given atom is controlled by a particular DOF.  Record
	/// that fact.
	virtual
	void
	add_atom(
		AtomID const & atom_id,
		DOF_ID const & dof_id
	);

	/// @brief Traverse the atom trees in preparation for minimization to tie together all the
	/// DOFs and the atoms they control.
	void
	setup( AtomTreeCollectionOP trees );

public:

	/// Accessors
	Size nactive_residues() const { return nactive_residues_; }
	Size active_residue( Size index ) const {debug_assert( index <= nactive_residues_ ); return active_residues_[ index ]; }

	/// @brief MinimizerMapBase class virtual accessor
	virtual kinematics::DomainMap const & domain_map() const { return domain_map_; }

	/// @brief Inline accessor
	inline kinematics::DomainMap const & dm() const { return domain_map_; }

	Size n_dof_nodes() const { return nactive_moving_atoms_total_; }

	/// @brief Initialize a multivec with the dofs reflected in the current residue(s)
	void starting_dofs( optimization::Multivec & dofs ) const;

	/// @brief Assign the chi values to the residue(s)
	void assign_dofs_to_mobile_residues( optimization::Multivec const & dofs );

	optimization::DOF_Node &
	dof_node( Size index );

	virtual
	conformation::Residue const &
	residue( Size seqpos ) const;

	optimization::DOF_Node const &
	dof_node_for_chi( Size resid, Size chiid ) const;

	id::TorsionID
	tor_for_dof( id::DOF_ID const & dofid ) const;

	kinematics::tree::Atom const &
	atom( AtomID const & atid ) const;

	void zero_atom_derivative_vectors();

	/// @brief propagate f1/f2's up from children to parents
	void link_torsion_vectors();

	void set_natoms_for_residue( Size resid, Size natoms );

	Size get_atom_index( id::AtomID const & atm ) {
	debug_assert (atm.rsd()>0 && atm.rsd()<=atoms_to_dofid_.size());
	debug_assert (atm.atomno()>0 && atm.atomno()<=atoms_to_dofid_[atm.rsd()].size());
		return atoms_to_dofid_[atm.rsd()][atm.atomno()];
	}

	id::AtomID const & get_atom( Size idx ) {
	debug_assert (idx>0 && idx<dofid_to_atoms_.size());
		return dofid_to_atoms_[idx];
	}

	optimization::MultifuncOP
	make_multifunc(
		pose::Pose & p,
		utility::vector1< conformation::ResidueCOP > const & bg_residues,
		scoring::ScoreFunction const & sfxn,
		scoring::MinimizationGraph & mingraph)
	{
		optimization::MultifuncOP retval( new CartSCMinMultifunc(p,bg_residues,sfxn,mingraph,*this) );
		return retval;
	}


protected:
	void reset_dof_nodes();

private:
	// cartesian dofs (split per-residue)
	utility::vector1< utility::vector1< id::AtomID > > moving_atoms_;
	utility::vector1< core::Size > nactive_moving_atoms_;
	core::Size nactive_moving_atoms_total_;

	// map atomIDs <-> dof indices
	utility::vector1< id::AtomID > dofid_to_atoms_;
	utility::vector1< utility::vector1< Size > > atoms_to_dofid_;

	// temporary storage for converting coords <-> multivec
	utility::vector1< core::Vector > residue_coord_workspace_;  // temporary space

	// even though we don't need _all_ of this data, use the same datatypes as atomtree variant
	AtomTreeCollectionOP atom_tree_collection_;
	utility::vector1< ResidueAtomTreeCollectionOP > atcs_for_residues_;

	kinematics::DomainMap domain_map_;
};

extern optimization::DOF_NodeOP dummy_nodeop;


} // namespace scmin
} // namespace pack
} // namespace core

#endif
