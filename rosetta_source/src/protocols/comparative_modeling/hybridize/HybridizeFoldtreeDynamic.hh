// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.hh
/// @author Yifan Song

#ifndef INCLUDED_protocols_comparative_modeling_hybridize_HybridizeFoldtreeDynamic_hh
#define INCLUDED_protocols_comparative_modeling_hybridize_HybridizeFoldtreeDynamic_hh

// Unit headers
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loops.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace hybridize {

class HybridizeFoldtreeDynamic {
public:
	HybridizeFoldtreeDynamic( );

	// initialize pose ... add VRT if needed
	void initialize( core::pose::Pose & pose, protocols::loops::Loops const & core_chunks );
    void reset( core::pose::Pose & pose );

private:
	// HELPER FUNCTIONS

	// stochastically selects an anchor position
	core::Size choose_anchor_position(const protocols::loops::Loop& chunk) const;

	// update pose fold tree with new chunks .. assume vrt already added
	void update(core::pose::Pose & pose);

	// set the core chunks ... where anchors must lie and cuts are forbidden
	void set_core_chunks(const protocols::loops::Loops & chunks);

	// cutpoints selection logic
	utility::vector1 < core::Size > decide_cuts(core::Size n_residues);

	// anchor selection logic
	void choose_anchors();

	void make_complete_chunks(
		utility::vector1 < core::Size > cut_positions,
		core::Size n_residues);

private:
	// DATA

	// segment the entire pose, defining cutpoints
	protocols::loops::Loops chunks_last_; 
	protocols::loops::Loops chunks_; 

	// segment pose into "core regions" where anchors may lie
	protocols::loops::Loops core_chunks_last_;
	protocols::loops::Loops core_chunks_;

	utility::vector1 < core::Size > anchor_positions_last_;
	utility::vector1 < core::Size > anchor_positions_;

	/// index of the virtual residue we added to the pose in set_up()
	int virtual_res_;
	core::Size num_nonvirt_residues_;

	// backup original info
	core::kinematics::FoldTree saved_ft_;
	core::Size saved_n_residue_;
};

}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols

#endif  // PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEDYNAMIC_HH_
