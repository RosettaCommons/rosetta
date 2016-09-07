// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/MinimizerMap.hh
/// @brief  Class for connecting DOFs in the atom tree to DOFs optimized by the AtomTreeMinimizer
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_MinimizerMap_hh
#define INCLUDED_core_optimization_MinimizerMap_hh

// Unit headers
#include <core/optimization/MinimizerMap.fwd.hh>

// Package headers
#include <core/optimization/DOF_Node.hh>
//#include <core/optimization/AtomNode.hh>
#include <core/optimization/types.hh>


// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/scoring/DerivVectorPair.hh>

// Rosetta headers
// #include <util_basic.hh>
// #include <jump_classes.hh>
// #include <core/kinematics/Stub.hh>
// #include <id.hh>

// // Numeric headers
// #include <numeric/all.fwd.hh>
// #include <numeric/conversions.hh>
// #include <numeric/xyzMatrix.hh>
// #include <numeric/xyzVector.hh>

// // ObjexxFCL headers
// #include <ObjexxFCL/FArray1D.hh>

// // Utility headers
// #include <utility/io/all.fwd.hh>

// // C++ headers
// #include <algorithm>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// //#include <iosfwd>
// #include <utility/assert.hh>
// #include <vector>
// #include <string>
// #include <map>
#include <list>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>


namespace core {
namespace optimization {


class MinimizerMap : public kinematics::MinimizerMapBase
{
public:
	typedef std::list< DOF_NodeOP > DOF_Nodes;
	typedef DOF_Nodes::iterator iterator;
	typedef DOF_Nodes::const_iterator const_iterator;

	typedef id::AtomID AtomID;
	typedef id::DOF_ID DOF_ID;
	typedef id::DOF_Type DOF_Type;
	typedef kinematics::MinimizerMapBase parent;

	typedef scoring::DerivVectorPair DerivVectorPair;

public:


	MinimizerMap():
		parent(),
		dof_node_pointer_()
	{}


	~MinimizerMap() override;


	void
	setup(
		pose::Pose & pose,
		kinematics::MoveMap const & move_map
	);



	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	) override;



	void
	add_atom(
		AtomID const & atom_id,
		DOF_ID const & dof_id
	) override;


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
	) const;



	kinematics::DomainMap const &
	domain_map() const override
	{
		return domain_map_;
	}

	/// get a dof_node by its dof_id
	DOF_NodeOP
	dof_node_from_id( DOF_ID const &id )
	{
		DOF_NodeOP node = nullptr;
		if ( id.valid() ) {
			node = dof_node_pointer_[ id ];
			if ( node == nullptr ) {
				std::cerr << "DOF_ID does not exist in map! torsion= " << id << std::endl;
				utility_exit();
			}
		}
		return node;
	}

	utility::vector1< DerivVectorPair > &
	atom_derivatives( Size resid ) {
		return atom_derivatives_[ resid ];
	}

private:

	/// deletes and clears dof_nodes_
	void
	clear_dof_nodes();

	void
	assign_rosetta_torsions( pose::Pose const & pose ); // part of setup

private:

	/// list of all the moving degrees of freedom
	DOF_Nodes dof_nodes_;

	/// pointer from DOF_ID to the corresponding DOF_NodeOP
	id::DOF_ID_Map< DOF_NodeOP > dof_node_pointer_;


	kinematics::DomainMap domain_map_;

	utility::vector1< utility::vector1< DerivVectorPair > > atom_derivatives_;

}; // MinimizerMap


} // namespace kinematics
} // namespace core


#endif
