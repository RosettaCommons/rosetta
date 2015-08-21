// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SCMinMinimizerMap.hh
/// @brief  Class for identifying the sidechain DOFs in the AtomTree which are free during
///         any particular call to the minimizer.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_SCMinMinimizerMap_hh
#define INCLUDED_core_pack_scmin_SCMinMinimizerMap_hh

// Unit headers
#include <core/pack/scmin/SCMinMinimizerMap.fwd.hh>

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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/optimization/types.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/optimization/Multifunc.fwd.hh>

//#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

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
class SCMinMinimizerMap : public kinematics::MinimizerMapBase
{
public:
	SCMinMinimizerMap() :
		focused_residue_( 0 ),
		nactive_residues_( 0 ),
		nonideal_(false)
	{}

	virtual
	~SCMinMinimizerMap() {}


	/// @brief the SCMinMinimizerMap has to know how many residues are in the pose; this allows
	/// it to do O(1) updates to its DomainMap -- this function costs O(N).
	virtual
	void set_total_residue( Size total_residue ) = 0;

	/// @brief Disable the minimization for all residues.  Ammortized O(1).
	virtual
	void clear_active_dofs() = 0;

	/// @brief Activate all the dofs for a particular residue.  Ammortized O(1).
	virtual
	void activate_residue_dofs( Size resindex ) = 0;

	/// @brief Invoked during the depth-first traversal through the AtomTree.  The AtomTree
	/// is indicating that a particular torsion is dependent on another torsion.  Record
	/// that fact.
	//fpd seperate implementations in cart and atomtree
	virtual
	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	) = 0;

	/// @brief Invoked during the depth-first traversal through the AtomTree; the atom
	/// tree is indicating that a given atom is controlled by a particular DOF.  Record
	/// that fact.
	//fpd seperate implementations in cart and atomtree
	virtual
	void
	add_atom(
		AtomID const & atom_id,
		DOF_ID const & dof_id
	) = 0;

	/// @brief Traverse the atom trees in preparation for minimization to tie together all the
	/// DOFs and the atoms they control.
	//fpd seperate implementations in cart and atomtree
	virtual
	void
	setup( AtomTreeCollectionOP trees ) = 0;

public:

	/// Accessors
	Size nactive_residues() const { return nactive_residues_; }
	Size active_residue( Size index ) const {debug_assert( index <= nactive_residues_ ); return active_residues_[ index ]; }

	/// @brief MinimizerMapBase class virtual accessor
	virtual kinematics::DomainMap const & domain_map() const { return domain_map_; }

	/// @brief Inline accessor
	inline kinematics::DomainMap const & dm() const { return domain_map_; }

	/// @brief Convenience lookup -- turns over the request to the AtomTreeCollection
	virtual
	conformation::Residue const &
	residue( Size seqpos ) const = 0;

	virtual
	Size n_dof_nodes() const = 0;

	/// @brief Initialize a multivec with the dofs reflected in the current residue(s)
	virtual
	void starting_dofs( optimization::Multivec & dofs ) const = 0;

	/// @brief Assign the chi values to the residue(s)
	virtual
	void assign_dofs_to_mobile_residues( optimization::Multivec const & dofs ) = 0;

	// fpd get dof_node by index
	virtual
	optimization::DOF_Node &
	dof_node( Size index ) = 0;

	virtual
	optimization::DOF_Node const &
	dof_node_for_chi( Size resid, Size chiid ) const = 0;

	virtual
	id::TorsionID
	tor_for_dof( id::DOF_ID const & dofid ) const = 0;

	virtual
	kinematics::tree::Atom const &
	atom( AtomID const & atid ) const = 0;

	virtual
	void zero_atom_derivative_vectors() = 0;

	/// @brief propagate f1/f2's up from children to parents
	virtual
	void link_torsion_vectors() = 0;

	virtual
	void set_natoms_for_residue( Size resid, Size natoms ) = 0;

	utility::vector1< scoring::DerivVectorPair > &
	atom_derivatives( Size resid ) {
		return atom_derivatives_[ resid ];
	}

	void
	set_nonideal( bool val_in ) {
		nonideal_ = val_in;
	}

	virtual
	optimization::MultifuncOP
	make_multifunc(
		pose::Pose & p,
		utility::vector1< conformation::ResidueCOP > const & bg_residues,
		scoring::ScoreFunction const & sfxn,
		scoring::MinimizationGraph & mingraph) = 0;

protected:
	virtual
	void reset_dof_nodes() = 0;

protected:
	/// each atom tree in the AtomTreeCollection will tell us that it represents residue 1.
	/// this variable tells us which residue is actually being represented.
	Size focused_residue_;

	Size nactive_residues_;

	utility::vector1< utility::vector1< scoring::DerivVectorPair > > atom_derivatives_;

	utility::vector1< Size > active_residue_index_for_res_;
	utility::vector1< Size > active_residues_;
	kinematics::DomainMap    domain_map_;

	// fpd for cartesian control bb movement instead
	bool nonideal_;
};


} // namespace scmin
} // namespace pack
} // namespace core

#endif
