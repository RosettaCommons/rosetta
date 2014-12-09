// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.hh
/// @brief  small class to setup an invrot tree in an existing pose
/// @author Florian Richter, flosopher@gmail.com, march 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_AlignPoseToInvrotTreeMover_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_AlignPoseToInvrotTreeMover_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// package headers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
//#include <utility/pointer/ReferenceCount.hh>

// C++ headers

//#include <utility/vector1.fwd.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief small mover that takes an invrot tree
class AlignPoseToInvrotTreeMover : public protocols::moves::Mover {

public:
	typedef core::Size Size;

	/// @brief the invrot tree and the seqpos object
	/// that are being passed in need to be initialized
	AlignPoseToInvrotTreeMover(
		InvrotTreeCOP invrot_tree,
		AllowedSeqposForGeomCstCOP seqpos
	);

	~AlignPoseToInvrotTreeMover();

	virtual
	std::string
	get_name() const;

	virtual
	void
	apply( core::pose::Pose & pose );

	void
	set_add_target_to_pose( bool const setting );

	void
	set_geomcst_for_superposition_from_enz_io(
		EnzConstraintIOCOP enzcst_io );

	/// @brief sets up a foldtree such that
	/// the anchor residue doesn't move,
	/// i.e. a backward edge from anchor to 1
	/// and a forward edge from anchor to seqpos
	/// also need to setup the jumps at anchor res
	void
	setup_foldtree_around_anchor_invrot(
		core::pose::Pose & pose,
		Size const anchor_seqpos,
		Size const first_target_seqpos ) const;

	/// @brief silly helper function to get the equivalent
	/// of a certain residue in a new residue type set
	/// note if the passed in residue is already of the
	/// required residue type set, nothing changes, just
	/// passes back input pointer
	core::conformation::ResidueCOP
	switch_residue_type_set(
		core::conformation::ResidueCOP residue,
		std::string const desired_restype_set_name
	) const;

private:

	bool add_target_to_pose_; //this variable decides whether the target residues are added or not
	InvrotTreeCOP invrot_tree_;
	AllowedSeqposForGeomCstCOP seqpos_;
	utility::vector1< InvrotCollectorCOP > all_invrots_;

	//the geomcsts for which rotamers the pose can be superimposed on
	utility::vector1< core::Size > geomcsts_for_superposition_;

};


}
}
}

#endif
