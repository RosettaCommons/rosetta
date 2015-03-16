// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hybridization/HybridizeFoldtreeDynamic.hh
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_HybridizeFoldtreeDynamic_hh
#define INCLUDED_protocols_hybridization_HybridizeFoldtreeDynamic_hh

// Unit headers
#include <protocols/hybridization/HybridizeFoldtreeDynamic.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loops.hh>

#include <utility/vector1.hh>

#include <set>

namespace protocols {
namespace hybridization {

class HybridizeFoldtreeDynamic {
public:
	HybridizeFoldtreeDynamic( );

	// initialize pose ... add VRT if needed
	void initialize(
				core::pose::Pose & pose,
				protocols::loops::Loops const & core_chunks,
				utility::vector1< std::pair< core::Size, core::Size > > const & strand_pairs,
				std::set<core::Size> const & strand_pair_library_positions
	);
	void reset( core::pose::Pose & pose );

	utility::vector1 < core::Size > get_anchors() { return anchor_positions_; }

private:
	// HELPER FUNCTIONS

	// stochastically selects an anchor position
	core::Size choose_anchor_position(const protocols::loops::Loop& chunk) const;

	// update pose fold tree with new chunks .. assume vrt already added
	void update(core::pose::Pose & pose);

	// set the core chunks ... where anchors must lie and cuts are forbidden
	void set_core_chunks(const protocols::loops::Loops & chunks);

	// fold tree
	//void jumps_and_cuts_from_pose( core::pose::Pose & pose, utility::vector1< std::pair<int, int > > & jumps, utility::vector1< int > & cuts);
	void jumps_and_cuts_from_foldtree( core::kinematics::FoldTree & foldtree, utility::vector1< std::pair<core::Size, core::Size > > & jumps, utility::vector1< core::Size > & cuts);

	// cutpoints selection logic
	utility::vector1 < core::Size > decide_cuts(core::pose::Pose & pose, core::Size n_residues);

	// anchor selection logic
	void choose_anchors();

	void make_complete_chunks(
		utility::vector1 < core::Size > cut_positions,
		core::Size n_residues);

	// straind pairings

	// gets core chunks that are strand paired to a specific chunk
	void get_pair_core_chunks( core::Size const chunk_index, utility::vector1<core::Size> & pair_chunks, utility::vector1<std::pair<core::Size, core::Size> > & pair_chunks_pairs );

	// gets a core chunk that covers a specific position
	void get_core_chunk_index_from_position( core::Size const position, core::Size & index );

	void add_overlapping_pair_chunks(
    core::Size const index,
    utility::vector1<core::Size> & cuts,
    utility::vector1<std::pair<core::Size, core::Size> > & jumps,
    std::set<core::Size> & rooted_chunk_indices );

	void remove_cut( core::Size const cut, utility::vector1<core::Size> & cuts );

private:
	// DATA

	// chunks cover the whole pose
	protocols::loops::Loops chunks_;

	// core chunks cover "rigid segments" (generally SS elts)
	protocols::loops::Loops core_chunks_;

	// anchor positions are within core chunks
	utility::vector1 < core::Size > anchor_positions_;

	/// index of the virtual residue we added to the pose in set_up()
	int virtual_res_;
	core::Size num_nonvirt_residues_;
	core::Size num_protein_residues_;

	// backup original info
	core::kinematics::FoldTree initial_asymm_foldtree_;
	core::kinematics::FoldTree saved_ft_;
	core::Size saved_n_residue_;

	// strand pairings
	utility::vector1< std::pair< core::Size, core::Size > > strand_pairs_; // res pos of strand pairs
	std::set<core::Size> strand_pair_library_positions_; // strand pair positions that are from the pairing library (not from a pdb template)
	std::set<core::Size> template_core_chunk_indices_; // core_chunks_ indices of chunks from templates (not from pairing library)
	std::set<core::Size> rooted_chunk_indices_; // chunk indices that have a jump to the root of the star fold tree

};

}  //  namespace hybridization
}  //  namespace protocols

#endif  // INCLUDED_protocols_hybridization_HybridizeFoldtreeDynamic_hh
