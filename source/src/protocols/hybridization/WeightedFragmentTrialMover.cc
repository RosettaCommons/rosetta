// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Fragment insertion and trial, each residue has a customized weight for the frequency of insertion
/// @author Yifan Song

#include <protocols/hybridization/WeightedFragmentTrialMover.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

#include <numeric/random/random.hh>
#include <utility/exit.hh>


namespace protocols {
namespace hybridization {

WeightedFragmentTrialMover::WeightedFragmentTrialMover(
	utility::vector1< core::fragment::FragSetOP > const frag_libs,
	utility::vector1< core::Real > const residue_weights,
	utility::vector1< core::Size > const anchor_residues,
	core::Size const nr_frags)
{
	moves::Mover::type( "WeightedFragmentTrialMover" );
	frag_libs_ = frag_libs;
	anchor_reses_ = anchor_residues;
	weighted_sampler_.resize(frag_libs.size());
	update_sampler_weights(residue_weights);
	nr_frags_ = nr_frags;
}

void WeightedFragmentTrialMover::update_sampler_weights( utility::vector1< core::Real > const residue_weights )
{
	total_frames_ = 0;
	for (Size i_frag_set = 1; i_frag_set<=frag_libs_.size(); ++i_frag_set) {
		utility::vector1< core::Real > frame_weights(frag_libs_[i_frag_set]->nr_frames(), 0.0);
		for (Size i_frame = 1; i_frame <= frag_libs_[i_frag_set]->nr_frames(); ++i_frame) {
			core::fragment::ConstFrameIterator frame_it = frag_libs_[i_frag_set]->begin(); // first frame of the fragment library
			advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			// disallow insertion across anchors
			bool cross_anchor = false;
			for (Size i_anchor = 1; i_anchor <= anchor_reses_.size() && !cross_anchor; ++i_anchor) {
				if (anchor_reses_[i_anchor]>=seqpos_start && anchor_reses_[i_anchor]<=seqpos_end) cross_anchor = true;
			}
			if (cross_anchor) continue;

			for (Size seqpos = seqpos_start; seqpos <= seqpos_end; ++seqpos) { // accumulate the weights of all residues in the fragment
				if (seqpos < 1 || seqpos > residue_weights.size()) {
					std::cerr << "seqpos = " << seqpos << " not in [1," << residue_weights.size() << "]" << std::endl;
					utility_exit_with_message("FATAL. Fragment library size doesn't match with the size of protein.");
				}
				frame_weights[i_frame] += residue_weights[seqpos];
			}
			if (frame_weights[i_frame] > 0.0) total_frames_++;
		}
		weighted_sampler_[i_frag_set].weights(frame_weights);
	}
}

void WeightedFragmentTrialMover::apply(core::pose::Pose & pose)
{
	// pick fragment set
	Size i_frag_set = numeric::random::rg().random_range(1, frag_libs_.size());
	// pick insertion position
	Size insert_pos = weighted_sampler_[i_frag_set].random_sample(numeric::random::rg());

	core::fragment::ConstFrameIterator frame_it = frag_libs_[i_frag_set]->begin();
	advance(frame_it, insert_pos-1);
	core::Size nr_frags = frame_it->nr_frags();
	if (nr_frags_ && nr_frags_ < nr_frags) nr_frags = nr_frags_;
	Size i_frag = numeric::random::rg().random_range(1, nr_frags);

	frame_it->apply( i_frag, pose );
}

std::string WeightedFragmentTrialMover::get_name() const
{
	return "WeightedFragmentTrialMover";
}

} // hybridization
} // protocols
