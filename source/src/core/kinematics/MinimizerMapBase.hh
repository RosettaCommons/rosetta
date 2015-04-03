// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/MinimizerMapBase.hh
/// @brief  Deriviative classes from this base class will implement three virtual functions
///         allowing the AtomTree to communicate with an outside Minimizer about how the
///         degrees of freedom in the AtomTree should be minimized.
/// @author Phil Bradley
/// @author Andrew Leaver-Fay -- refactoring phil's code a little.

#ifndef INCLUDED_core_kinematics_MinimizerMapBase_hh
#define INCLUDED_core_kinematics_MinimizerMapBase_hh

// Unit headers
#include <core/kinematics/MinimizerMapBase.fwd.hh>

// Package headers
//#include <core/optimization/DOF_Node.hh>
//#include <core/optimization/AtomNode.hh>
//#include <core/optimization/types.hh>


// Project headers
//#include <core/pose/Pose.fwd.hh>
//#include <core/id/AtomID_Map.hh>
//#include <core/id/DOF_ID_Map.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/DomainMap.fwd.hh>


namespace core {
namespace kinematics {


class MinimizerMapBase : public utility::pointer::ReferenceCount
{
public:

	typedef id::AtomID AtomID;
	typedef id::DOF_ID DOF_ID;
	typedef id::DOF_Type DOF_Type;

public:

	/// @brief default ctor; noop
	MinimizerMapBase();

	/// @brief dstor
	virtual ~MinimizerMapBase();


	/// @brief Allow the AtomTree to communicate to this class that a particular torsion
	/// (or angle or distance -- a particular DOF) belongs in the minimization task
	/// to inform this class the parent DOF for that torsion.
	virtual
	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	) = 0;

	/// @brief Allow the AtomTree to inform this class that a particular atom belongs
	/// in the derivative calculation for a certain DOF.  That certain DOF must have already
	/// been included in the minimization process through a prior invocation of the
	/// add_torsion method.
	virtual
	void
	add_atom(
		AtomID const & atom_id,
		DOF_ID const & dof_id
	) = 0;


	/*
	const_iterator
	begin() const
	{
		return dof_nodes_.begin();
	}


	const_iterator
	end() const
	{
		return dof_nodes_.end();
	}


	iterator
	begin()
	{
		return dof_nodes_.begin();
	}


	iterator
	end()
	{
		return dof_nodes_.end();
	}


	DOF_Nodes const &
	dof_nodes() const
	{
		return dof_nodes_;
	}


	DOF_Nodes &
	dof_nodes()
	{
		return dof_nodes_;
	}

	/// this will only work if DOF_nodes are sorted!
	void
	link_torsion_vectors();


	void
	zero_torsion_vectors();

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


	inline
	int
	nangles() const
	{
		return dof_nodes_.size();
	}


	void
	reset_jump_rb_deltas(
		pose::Pose & pose,
		Multivec & dofs
	) const;

	// this is used in pack/unpack_phipsi and deriv calculation
	Real
	torsion_scale_factor(
		DOF_Node const & tor
	) const; */


	virtual
	kinematics::DomainMap const &
	domain_map() const = 0;

	/*
	/// get a dof_node by its dof_id
	DOF_Node*
	dof_node_from_id( DOF_ID const &id )
	{
		DOF_Node* node = 0;
		if ( id.valid() ) {
			node = dof_node_pointer_[ id ];
			if ( node == 0 ) {
				std::cerr << "DOF_ID does not exist in map! torsion= " << id << std::endl;
				utility_exit();
			}
		}
		return node;
	}

private:

	/// deletes and clears dof_nodes_
	void
	clear_dof_nodes();

	/// list of all the moving degrees of freedom
	DOF_Nodes dof_nodes_;

	/// pointer from DOF_ID to the corresponding DOF_Node*
	id::DOF_ID_Map< DOF_Node* > dof_node_pointer_;


	kinematics::DomainMap domain_map_;

	void
	assign_rosetta_torsions( pose::Pose const & pose ); // part of setup

*/
}; // MinimizerMapBase


} // namespace kinematics
} // namespace core


#endif
