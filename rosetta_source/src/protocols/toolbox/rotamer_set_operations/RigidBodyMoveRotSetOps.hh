// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/graph/Graph.fwd.hh>
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
/// body movements in the packer. it holds a list of alternative rigid body
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
/// In the base class, the alternative positions can only be set from the outside.
/// there are many other conceivable ways of generating alternative confs though
/// (docking, systematically, etc). Ideally, if one wanted to generate them in this
/// way, one would need to write a new class derived from this base class that
/// only deals with generating the alternative rigid body conf, and this base class
/// handles the packer specific stuff (superposition/neighbor radius)
class RigidBodyMoveRSO : public core::pack::rotamer_set::RotamerSetOperation
{
public:
	typedef core::pack::rotamer_set::RotamerSetOperation parent;
	typedef core::Real Real;
	typedef core::Size Size;

	RigidBodyMoveRSO( core::Size seqpos );
	RigidBodyMoveRSO( RigidBodyMoveRSO const & other );
	~RigidBodyMoveRSO();

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::pack::task::PackerTask const & ptask,
		core::graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);

	virtual
	Real
	increase_packer_residue_radius(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP //the_task
	) const;

	void
	set_rigid_body_confs(
		utility::vector1< core::conformation::ResidueCOP > const & rigid_body_confs
	);

	/// @brief returns the largest observed distance between the nbr atom
	/// in the target res and the nbr atom in any of the saved rb confs
	Real
	determine_largest_nbr_atom_distance(
		core::conformation::Residue const & target_res
	) const;


protected:


private:

	core::Size seqpos_;
	utility::vector1< core::conformation::ResidueCOP > rigid_body_confs_;
	core::pack::rotamer_set::BumpSelector bump_selector_;
};

} //namespace protocols
} //namespace toolbox
} //namespace rotamer_set_operations

#endif

