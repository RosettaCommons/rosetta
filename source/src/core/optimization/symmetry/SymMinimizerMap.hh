// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/SymMinimizerMap.hh
/// @brief  MinimizerMap for symmetric minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_optimization_symmetry_SymMinimizerMap_hh
#define INCLUDED_core_optimization_symmetry_SymMinimizerMap_hh

// Unit headers
#include <core/optimization/symmetry/SymMinimizerMap.fwd.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/DOF_Node.fwd.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <list>

#include <core/scoring/DerivVectorPair.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>

#ifdef WIN32
	#include <core/optimization/DOF_Node.hh>
#endif

namespace core {
namespace optimization {
namespace symmetry {

/// @brief Atom tree multifunction class
class SymMinimizerMap : public kinematics::MinimizerMapBase {
public:
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;
	typedef std::list< DOF_NodeOP >                 DOF_Nodes;
	typedef DOF_Nodes::const_iterator               const_iterator;
	//typedef DOF_Nodes::iterator                     iterator;

public: // Creation

	// c-tor
	SymMinimizerMap(
		pose::Pose const & pose, // must have been scored before this setup; energy map must be up-to-date
		kinematics::MoveMap const & mm, // does not have to be "symmetric" -- it will be symmetrized
		SymmetryInfoCOP symm_info,
		bool const new_sym_min = false
	);

	/// @brief Destructor
	virtual
	~SymMinimizerMap();

	/// @brief The atom tree will report that a new torsion has been identified as free in the traversal of the atom tree.
	/// If this is an independent torsion, then the SymMinimizerMap will add a new DOF_Node, but otherwise, will
	/// ignore the DOF.  The atom tree will traverse through dependent torsions in addition to independent torsions, and
	/// it's the job of the SymMinimizerMap to weed out the dependent torsions.
	virtual
	void
	add_torsion(
		DOF_ID const & new_torsion,
		DOF_ID const & parent
	);

	/// @brief Add an atom to the list of atoms controlled by a given DOF.  The SymMinimzierMap
	/// will figure out, first, if the dof_id is a dependent or independent dof.  If it's a dependent
	/// DOF, then it will figure out if the given atom has any interactions with an independent residue.
	/// If not, then the atom is ignored.  If it does, then the SymMinimizerMap will figure out
	/// what independent DOF the given dependent DOF is a a clone of, and add this atom as being controlled
	/// by that dependent DOF.
	virtual
	void
	add_atom(
		AtomID const & atom_id,
		DOF_ID const & dof_id
	);

	virtual
	kinematics::DomainMap const &
	domain_map() const;

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

	DOF_NodeOP
	dof_node_from_id( DOF_ID const &id ) const;

	Size nangles() const { return n_independent_dof_nodes_; }

	void zero_torsion_vectors();

	void link_torsion_vectors();

	Real
	torsion_scale_factor(
		DOF_Node const & dof_node
	) const;

	void
	reset_jump_rb_deltas(
		pose::Pose & pose,
		Multivec & dofs
	) const;

public:
	/// Allow read/write access to the DOF_Nodes themselves, but do not allow anyone to change the
	/// DOF_Nodes list.  Elements cannot be dropped from the list, nor should the list be clearable.
	/// Of course: if you have a const iterator to a list element containing a pointer, then it is
	/// entirely possible to perform non-const operations on the thing being pointed at.
	/// HOWEVER, neither the pointer nor the list element can be changed.

	/// @brief begin iterator for the independent dofs
	const_iterator
	begin() const
	{
		return dof_nodes_.begin();

	}
	/// @brief End iterator for the independent dofs
	const_iterator
	end() const
	{
		return dof_nodes_.end();
	}

	const_iterator
	dependent_begin() const
	{
		return dependent_dof_nodes_.begin();
	}

	///
	const_iterator
	dependent_end() const
	{
		return dependent_dof_nodes_.end();
	}



	///
	DOF_Nodes const &
	dof_nodes() const
	{
		return dof_nodes_;
	}

	/// @brief Retrieve the per-atom derivatives that are accumulated in to
	utility::vector1< scoring::DerivVectorPair > &
	atom_derivatives( Size resid ) {
		return atom_derivatives_[ resid ];
	}


	bool
	new_sym_min() const { return new_sym_min_; }

	/// @brief Convert a cloned dof into its equivalent in the asymmetric unit
	id::DOF_ID asymmetric_dof( DOF_ID const & cloned_dof ) const;

private:

	/// @brief Non-virtual method -- not invoked directly by the atom tree.
	void
	add_new_dof_node(
		DOF_ID const & new_torsion,
		DOF_ID const & parent,
		bool dependent
	);

	void assign_rosetta_torsions();


private:
	pose::Pose const & pose_; // needed in the add_torsion, and add_atom callback functions
	SymmetryInfoCOP symm_info_;
	utility::vector1< bool > res_interacts_with_asymmetric_unit_;

	DOF_Nodes dof_nodes_;
	DOF_Nodes dependent_dof_nodes_;
	//Size n_dof_nodes_; // this was not used
	Size n_independent_dof_nodes_;

	/// pointer from DOF_ID to the corresponding DOF_NodeOP
	id::DOF_ID_Map< DOF_NodeOP > dof_node_pointer_;

	kinematics::DomainMap domain_map_;

	utility::vector1< utility::vector1< scoring::DerivVectorPair > > atom_derivatives_;

	/// adding this guy so we can tell more accurately which dof's are dependent/independent
	id::DOF_ID_Map< id::TorsionID > dof_id2torsion_id_;

	///
	bool new_sym_min_;

};

} // symmetry
} // namespace optimization
} // namespace core

#endif
