// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Align a random jump to template
/// @details
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_InsertChunkMover_hh
#define INCLUDED_protocols_hybridization_InsertChunkMover_hh

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

#include <set>

namespace protocols {
namespace hybridization {

class InsertChunkMover: public protocols::moves::Mover
{

public:

	InsertChunkMover();
	~InsertChunkMover() override;

	void init();

	void set_bb_xyz_aligned(core::pose::Pose & pose);

	bool success() {
		return success_;
	}

	void set_template(core::pose::PoseCOP template_pose, core::Size template_id,
		std::map <core::Size, core::Size> const & sequence_alignment);

	void set_aligned_chunk(core::pose::Pose const & pose, Size const jump_number, bool anchor_insert_only_in);

	bool get_local_sequence_mapping(core::pose::Pose & pose,
		int registry_shift = 0,
		Size MAX_TRIAL = 100 );

	void check_overlap(core::pose::Pose & pose);

	void set_registry_shift(int registry_shift);

	Size trial_counter(Size ires);

	void apply(core::pose::Pose & pose) override;

	std::string get_name() const override;

private:
	core::pose::PoseCOP template_pose_;
	core::Size template_id_;
	std::map <core::Size, core::Size> sequence_alignment_;
	int registry_shift_;
	bool anchor_insert_only_;

	Size jump_number_; // the jump to be realigned
	Size seqpos_start_; // start and end seqpose of the chunk, downstream of the jump
	Size seqpos_stop_;
	Size seqpos_aligned_start_; // start and end seqpose of the aligned piece
	Size seqpos_aligned_stop_;

	bool success_;

	// parameters of the protocol
	bool align_to_ss_only_; // only use the secondary structure portion to align to the template
	bool copy_ss_torsion_only_; // only copy the secondary structure information from the template

	char secstruct_;

	std::map <core::Size, core::Size> sequence_alignment_local_; // with registry shift of the aligned chunk
	core::id::AtomID_Map< core::id::AtomID > atom_map_; // atom map for superposition
	utility::vector1 <Size> align_trial_counter_;

};

} // hybridization
} // protocols

#endif
