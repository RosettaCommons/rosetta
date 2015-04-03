// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/CartesianMinimizerMap.hh
/// @brief  Class for connecting DOFs in the atom tree to DOFs optimized by the AtomTreeMinimizer
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_CartesianMinimizerMap_hh
#define INCLUDED_core_optimization_CartesianMinimizerMap_hh

// Unit headers
#include <core/optimization/CartesianMinimizerMap.fwd.hh>

// Package headers
#include <core/optimization/types.hh>

#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/AtomID.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/DerivVectorPair.hh>


#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>


namespace core {
namespace optimization {

class CartesianMinimizerMap : public kinematics::MinimizerMapBase
{
public:


	CartesianMinimizerMap()
	{}


	~CartesianMinimizerMap();


	virtual
	kinematics::DomainMap const &
	domain_map() const
	{
		return domain_map_;
	}

	virtual
	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	);

	virtual
	void
	add_atom(
		AtomID const & AtomID,
		DOF_ID const & dof_id
	);


	void
	setup(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map
	);

	/// clears old data+dimensions dof_node_pointer using size data from the pose
	void
	reset( pose::Pose const & pose );


	void
	copy_dofs_from_pose(
		pose::Pose const & pose,
		Multivec & dofs
	) const;


	void
	copy_dofs_to_pose(
		pose::Pose & pose,
		Multivec const & dofs
	) const;


	inline int
	ndofs() const	{ return 3*moving_atoms_.size(); }

	inline int
	ntorsions() const	{ return moving_dofids_.size(); }

	inline id::DOF_ID
	get_dof_id(Size n) const { return moving_dofids_[n]; }

	inline id::TorsionID
	get_TorsionID(Size n) const { return moving_torsionids_[n]; }

	inline int
	natoms() const	{ return moving_atoms_.size(); }

	inline id::AtomID
	get_atom(Size n) const { return moving_atoms_[n]; }

	inline Size
	get_atom_index(id::AtomID x) const { return atom_indices_[x]; }

	utility::vector1< core::scoring::DerivVectorPair > &
	atom_derivatives( Size resid ) {
		return atom_derivatives_[ resid ];
	}

	void
	zero_stored_derivs();

private:

	// maps dofIDs->torsionIDs
	// removes dofIDs with no relevant atoms
	void
	assign_rosetta_torsions_and_trim( pose::Pose const & pose ); // part of setup


private:

	// map moving atoms->index
	id::AtomID_Map< core::Size > atom_indices_;

	/// list of all the moving atoms
	utility::vector1<id::AtomID> moving_atoms_;


	utility::vector1< utility::vector1< core::scoring::DerivVectorPair > > atom_derivatives_;

	/// list of all moving torsions: dof ids and torsion ids
	/// we don't need all the info from dof_node, just this
	utility::vector1<id::DOF_ID> moving_dofids_;
	utility::vector1<id::TorsionID> moving_torsionids_;


	kinematics::DomainMap domain_map_;

}; // CartesianMinimizerMap


} // namespace kinematics
} // namespace core


#endif
