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

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.fwd.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>


namespace core {
namespace pose {
namespace rna {

class RNA_IdealCoord {
public:

	RNA_IdealCoord();
	~RNA_IdealCoord();

	//Apply ideal coords to one residue. Keep the backbone torsion values by default
	void apply( Pose & pose, Size const seqpos, bool const is_north = true, bool const keep_backbone_torsion = true ) const;

	//Apply ideal coords to whole pose.
	//pucker_conformations: 0 for skipping, 1 for North, 2 for South
	void apply( Pose & pose, utility::vector1 < Size > const & pucker_conformations, bool const keep_backbone_torsion = true ) const;

private:
	void init();
	bool is_torsion_exists(Pose const & pose, id::TorsionID const & torsion_id) const;
	utility::vector1 < Pose > ref_pose_list_;
	std::string const path_;
};

}
}
}

#endif
