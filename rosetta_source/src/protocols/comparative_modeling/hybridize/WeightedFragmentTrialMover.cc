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

#include <protocols/comparative_modeling/hybridize/WeightedFragmentTrialMover.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG(8401848);

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

WeightedFragmentTrialMover::WeightedFragmentTrialMover(
	utility::vector1< core::fragment::FragSetOP > const frag_libs,
	utility::vector1< core::Real > const residue_weights,
	utility::vector1< core::Size > const anchor_residues )
{
	frag_libs_ = frag_libs;
	anchor_reses_ = anchor_residues;
	weighted_sampler_.resize(frag_libs.size());
	update_sampler_weights(residue_weights);
}

void WeightedFragmentTrialMover::update_sampler_weights( utility::vector1< core::Real > const residue_weights )
{
	for (Size i_frag_set = 1; i_frag_set<=frag_libs_.size(); ++i_frag_set) {
		utility::vector1< core::Real > frame_weights(frag_libs_[i_frag_set]->nr_frames(), 0.0);
		for (Size i_frame = 1; i_frame <= frag_libs_[i_frag_set]->nr_frames(); ++i_frame) {
			core::fragment::FrameIterator frame_it = frag_libs_[i_frag_set]->begin(); // first frame of the fragment library
			advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			// disallow insertion across anchors
			bool cross_anchor = false;
			for (Size i_anchor = 1; i_anchor <= anchor_reses_.size() && !cross_anchor; ++i_anchor) {
				if (i_anchor>=seqpos_start && i_anchor<=seqpos_end) cross_anchor = true;
			}
			if (cross_anchor) continue;

			for (Size seqpos = seqpos_start; seqpos <= seqpos_end; ++seqpos) { // accumulate the weights of all residues in the fragment
				frame_weights[i_frame] += residue_weights[seqpos];
			}
		}
		weighted_sampler_[i_frag_set].weights(frame_weights);
	}
}

void WeightedFragmentTrialMover::apply(core::pose::Pose & pose)
{
	// pick fragment set
	Size i_frag_set = RG.random_range(1, frag_libs_.size());
	// pick insertion position
	Size insert_pos = weighted_sampler_[i_frag_set].random_sample(RG);

	core::fragment::FrameIterator frame_it = frag_libs_[i_frag_set]->begin();
	advance(frame_it, insert_pos-1);
	Size i_frag = RG.random_range(1, frame_it->nr_frags());
	
	frame_it->apply( i_frag, pose );
}
	
std::string WeightedFragmentTrialMover::get_name() const
{
	return "WeightedFragmentTrialMover";
}

} // hybridize 
} // comparative_modeling 
} // protocols
