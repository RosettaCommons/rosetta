// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/AtomTreeCollection.hh
/// @brief  Classes for holding sets of AtomTrees used during variuos packing+minimizing schemes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_AtomTreeCollection_hh
#define INCLUDED_core_pack_scmin_AtomTreeCollection_hh

// Unit headers
#include <core/pack/scmin/AtomTreeCollection.fwd.hh>

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/DOF_ID.hh>


#ifdef WIN32
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/Residue.hh>
#endif

//#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzVector.hh>

namespace core {
namespace pack {
namespace scmin {

/// @brief Class to save the state for a ResidueAtomTreeCollection
class ResidueAtomTreeCollectionMomento : public utility::pointer::ReferenceCount
{
public:
	ResidueAtomTreeCollectionMomento();
	virtual ~ResidueAtomTreeCollectionMomento();

	ResidueAtomTreeCollectionMomento( ResidueAtomTreeCollectionMomento const & );
	ResidueAtomTreeCollectionMomento & operator = ( ResidueAtomTreeCollectionMomento const & );

	void set_restype_index( Size setting );
	void copy_coords( conformation::Residue const & );

	Size restype_index() const { return restype_index_; }

	Vector const &
	coord( Size ind ) const
	{
		debug_assert( ind <= natoms_ );
		return coords_[ ind ];
	}

private:
	Size restype_index_;
	Size natoms_;
	utility::vector1< Vector > coords_;
};

/// @brief The conformation::Residues and kinematics::AtomTrees for a single
/// collection of rotamers (e.g. a RotamerSet).  Each chemical::ResidueType gets its own residue/atomtree pair.
/// A particular AtomTree/Residue pair can be set as "active" and manipulated
/// by setter and getters for either the coordinates of the Residues or the
/// chi dihedrals of the AtomTree.
class ResidueAtomTreeCollection : public utility::pointer::ReferenceCount
{
public:
	ResidueAtomTreeCollection(
		task::ResidueLevelTask const & rltask,
		conformation::Conformation const & conformation,
		conformation::Residue const & orig_res
	);

	/// @brief Constructor for a RotamerSet that has already been built.
	ResidueAtomTreeCollection(
		rotamer_set::RotamerSet const & rset,
		Size resid
	);

	virtual ~ResidueAtomTreeCollection();

	Size active_restype_index() const { return active_restype_; }
	void set_active_restype_index( Size restype_index );

	/// @brief The responsibility for making sure that the active residue and the active atomtree
	/// are in synch is offloaded to an external class so that the calls to "active_residue()" and
	/// "active_atom_tree()" are as fast as possible (and bitwise const for future multithreaded use).
	/// After a round of set_chi() calls, the user for this class must update the residue coordinates.
	void update_residue();

	/// @brief See comments for update_residue().  After a call to "set_rescoords", the user must
	/// call update_atomtree() to make sure the atomtree and the residue are in synch.
	void update_atom_tree();

	chemical::ResidueType const & active_restype() const;

	/// @brief Inline accessor for the active residue.
	/// asserts residue_uptodate_ -- make sure that update_residue is called first
	inline
	conformation::Residue const &  active_residue() const
	{
		debug_assert( residue_uptodate_ );
		return *active_residue_;
	}

	conformation::ResidueCOP active_residue_cop() const;

	/// @brief Inline accessor for the active atom tree.
	/// asserts atom_tree_uptodate_ -- make sure that update_atom_tree is called first
	kinematics::AtomTree const & active_atom_tree() const
	{
		debug_assert( atom_tree_uptodate_ );
		return *active_atom_tree_;
	}

	///fpd  get dof value by DOF_ID
	///     NOTE input resid must be 1
	core::Real dof( core::id::DOF_ID const &dofid );

	/// @brief Assigns the chi dihedral for the active restype.  Must be followed by a call to
	/// update_residue() before the next call to active_residue()
	void set_chi( Size chi_index, Real value );

	///fpd bondlength analog to set_chi
	///    like set_chi, assumes changes propagate to atomtree
	///    keyed off of chi#, so we only allow distances corresponding to chi angles to refine
	///    distance corresponds to the distance between atoms 3 and 4 defining the chi
	///    chino==0 ==> CA-CB distance, which allows us to refine ALA CB position for example
	void set_d( Size chi_index, Real value );

	///fpd bondangle analog to set_chi
	///    same idea as set_d
	void set_theta( Size chi_index, Real value );

	/// @brief Assigns the coordinates for a residue.  Must be followed by a call to
	/// update_atom_tree() before the next cal to active_atom_tree().
	void set_rescoords( conformation::Residue const & res );
	void set_rescoords( utility::vector1< Vector > const & coords );
	void set_rescoords( utility::vector1< id::AtomID > const & atms, utility::vector1< Vector > const & coords );

	void save_momento( ResidueAtomTreeCollectionMomento & momento ) const;
	void update_from_momento( ResidueAtomTreeCollectionMomento const & momento );

private:

	//uncopyable -- unimplemented
	ResidueAtomTreeCollection( ResidueAtomTreeCollection const & );
	ResidueAtomTreeCollection & operator = ( ResidueAtomTreeCollection const & );


private:

	Size active_restype_;
	kinematics::AtomTreeOP active_atom_tree_;
	conformation::ResidueOP active_residue_;

	bool residue_uptodate_;
	bool atom_tree_uptodate_;

	utility::vector1< kinematics::AtomTreeOP > atom_tree_representatives_;
	utility::vector1< conformation::ResidueOP  > residue_representatives_;

};

/// @brief A collection of ResidueAtomTreeCollection objects for an entire design task.
class AtomTreeCollection : public utility::pointer::ReferenceCount
{
public:
	/// @brief construct a forest
	AtomTreeCollection(
		pose::Pose const & pose,
		rotamer_set::RotamerSets const &
	);

	/// @brief construct a forest -- only use when no RotamerSets is constructed
	AtomTreeCollection(
		pose::Pose const & pose,
		task::PackerTask const &
	);

	/// @brief construct a single ResidueAtomTreeCollection from an existing RotamerSet
	AtomTreeCollection(
		pose::Pose const & pose,
		rotamer_set::RotamerSet const & rset,
		Size resid
	);

	/// @brief construct a single ResidueAtomTreeCollection -- only use when no RotamerSet is constructed
	AtomTreeCollection(
		pose::Pose const & pose,
		task::ResidueLevelTask const & rltask,
		Size resid
	);

	virtual ~AtomTreeCollection();

	ResidueAtomTreeCollection & moltenres_atomtree_collection( Size moltenresid );
	ResidueAtomTreeCollection & residue_atomtree_collection( Size resid );
	ResidueAtomTreeCollectionOP residue_atomtree_collection_op( Size resid );

	//conformation::Residue const & residue( Size seqpos );

private:

	//rotamer_set::RotamerSetsCOP rotsets_;
	utility::vector1< Size > resid_2_moltenresid_;
	utility::vector1< Size > moltenresid_2_resid_;
	utility::vector1< ResidueAtomTreeCollectionOP > res_collections_;

};

} // namespace scmin
} // namespace pack
} // namespace core

#endif
