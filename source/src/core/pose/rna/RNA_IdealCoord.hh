// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/RNA_IdealCoord.hh
/// @brief Apply ideal RNA geometry to a residue or a pose
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_pose_rna_RNA_IdealCoord_HH
#define INCLUDED_core_pose_rna_RNA_IdealCoord_HH

// Unit headers
#include <core/pose/rna/RNA_IdealCoord.fwd.hh>
#include <core/id/DOF_ID_Map.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#ifdef WIN32
	#include <core/pose/MiniPose.hh>
#else
	#include <core/pose/MiniPose.fwd.hh>
#endif
// Project headers
#include <core/chemical/rna/RNA_Util.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

using core::chemical::rna::PuckerState;

namespace core {
namespace pose {
namespace rna {

class RNA_IdealCoord : public utility::pointer::ReferenceCount {
public:

	RNA_IdealCoord();
	~RNA_IdealCoord();

	//Apply ideal coords to one residue. Keep the backbone torsion values by default
	//std::map < id::DOF_ID , Real > apply_and_return( Pose & pose, Size const seqpos, Size pucker, bool const keep_backbone_torsion = true ) const;

	void apply( Pose & pose, Size const seqpos, PuckerState pucker, bool const keep_backbone_torsion = true ) const;

	//Apply ideal coords to whole pose.
	//pucker_conformations: 0 for keeping pucker, 1 for North, 2 for South
	void apply( Pose & pose, utility::vector1 < PuckerState > const & puckers, bool const keep_backbone_torsion = true ) const;

	//Apply ideal coords to whole pose. Keep all pucker state
	void apply( Pose & pose, bool const keep_backbone_torsion = true ) const;

	// Apply ideal coord to puckers only, assuming ideal coords of A
	void apply_pucker( Pose & pose, Size const seqpos, PuckerState pucker, bool const keep_backbone_torsion = true ) const;

	bool is_torsion_exists(Pose const & pose, id::TorsionID const & torsion_id) const;

private:
	void init();
	void apply_coords( Pose & pose, Size const seqpos, Size const res_class, bool const ignore_base, bool const keep_backbone_torsion ) const;
	//	utility::vector1 < PoseOP > ref_pose_list_;
	utility::vector1 < MiniPoseOP > ref_mini_pose_list_;
	std::string const path_;
	Real delta_cutoff_;
};

}
}
}

#endif
