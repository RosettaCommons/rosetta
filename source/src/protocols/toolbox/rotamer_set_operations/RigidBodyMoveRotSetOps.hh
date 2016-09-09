// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/RigidBodyMoveRotSetOps.hh
/// @brief  classes for rigid body movement during rotamer packing
/// @author Florian Richter, floric@u.washington.edu, sep 2009

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_RigidBodyMoveRotSetOps_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_RigidBodyMoveRotSetOps_hh

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/RigidBodyMoveRotSetOps.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

//Project headers
#include <core/conformation/Residue.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/BumpSelector.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

/// @details base class that handles single residue (usually ligand) rigid
/// body movements in the packer. subclass generates a list of alternative rigid body
/// placements. when the alter_rotamer_set function is called, this class
/// superimposes all rotamers in the initial RotamerSet onto every alternative
/// rigid body conformation in the member data list. these superimposed
/// conformations then get added to the rotamer set.
/// this class also handles the increased packer radius necessitated by the
/// rigid body movements. when the increase_packer_residue_radius function
/// is called, the biggest distance between the residue in the initial pose
/// and any of the alternate conformations is returned.
///
/// Generating the alternative rigid body positions:
///
///   Subclasses must implement get_rigid_body_confs.
///
///   Subclasses may implement increase_packer_residue_radius to
/// optimize this calculation. Default implementation calls get_rigid_body_confs
/// and calculates the maximum nbr_atom delta.
///
class RigidBodyMoveBaseRSO : public core::pack::rotamer_set::RotamerSetOperation
{
	typedef core::pack::rotamer_set::RotamerSetOperation parent;

protected:
	RigidBodyMoveBaseRSO() :
		parent()
	{}

	RigidBodyMoveBaseRSO(const RigidBodyMoveBaseRSO& other) :
		parent(other),
		bump_selector_(other.bump_selector_)
	{}


public:
	/// @brief Adds additional rotamers at each rb conf.
	///
	/// @details
	/// fairly simple: iterate over the rotamers in the rotamer_set, superimpose
	/// each of them onto all the internally stored rigid body confs, and then
	/// add the newly generated rotamers to the rotamer_set
	///
	/// Recalculates chi-angle expansion via call to RotamerSet::extra_chi_samples
	/// and then generate candidate rotamers from each rb conf via
	/// SingleResidueRotamerLibrary::fill_rotamer_vector if rotamer library is
	/// available for the residue type, else just adds the alternate confs.
	///
	/// Does not have full safety checks to make sure all the rotamers in the set are
	/// of the same residue type as the internally stored ones, i.e. that no user sets this
	/// position to designing in the task
	virtual
	void
	alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::pack::task::PackerTask const & ptask,
		utility::graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);

	/// @brief returns the largest possible change in nbr atoms location
	///
	/// Default implementation calls determine_larget_nbr_atom_distance with
	/// results of get_rigid_body_confs. Override this function if calculation
	/// can be optimized without call call to get_rigid_body_confs.
	using core::pack::rotamer_set::RotamerSetOperation::increase_packer_residue_radius;
	virtual
	core::Real
	increase_packer_residue_radius(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP, //the_task
		core::Size residue_index
	);

	/// @brief returns candidate alternate RB conformations
	virtual
	utility::vector1< core::conformation::ResidueCOP >
	get_rigid_body_confs(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & ptask,
		core::Size residue_index) = 0;

	/// @brief returns the largest observed distance between the nbr atom
	/// in the target res and the nbr atom in any of the candidate rb confs
	static
	core::Real
	determine_largest_nbr_atom_distance(
		core::conformation::Residue const & target_res,
		utility::vector1< core::conformation::ResidueCOP > alternate_confs);

private:
	core::pack::rotamer_set::BumpSelector bump_selector_;
};

/// @details
/// Basic implementation of alternate rb conf set operation.
/// The alternative positions are set externally before packing.
class RigidBodyMoveRSO : public RigidBodyMoveBaseRSO
{
public:
	typedef RigidBodyMoveBaseRSO parent;

	// 'seqpos' could be removed, it can be derived from the input
	// rotamer set. It is left in as a sanity check for this class.
	RigidBodyMoveRSO( core::Size seqpos );
	RigidBodyMoveRSO( RigidBodyMoveRSO const & other );

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	utility::vector1< core::conformation::ResidueCOP >
	get_rigid_body_confs(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & ptask,
		core::Size residue_index);

	void
	set_rigid_body_confs(
		utility::vector1< core::conformation::ResidueCOP > const & rigid_body_confs
	);

private:

	core::Size seqpos_;
	utility::vector1< core::conformation::ResidueCOP > rigid_body_confs_;
};

} //namespace protocols
} //namespace toolbox
} //namespace rotamer_set_operations

#endif

