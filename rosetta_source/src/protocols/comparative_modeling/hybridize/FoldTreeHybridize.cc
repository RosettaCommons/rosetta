// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef apps_pilot_yfsong_FoldTreeHybridize_HH
#define apps_pilot_yfsong_FoldTreeHybridize_HH

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <apps/pilot/yfsong/FoldTreeHybridize.fwd.hh>

#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>


namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::loops;
using namespace protocols::nonlocal;
using namespace numeric::model_quality;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;

	
public:
FoldTreeHybridize::FoldTreeHybridize(numeric::random::RandomGenerator & RG,
							 utility::vector1 < core::pose::PoseOP > const & template_poses,
							 Loops ss_chunks_pose,
							 AlignOption align_option = all_chunks,
							 Size max_registry_shift = 0) :
RG_(RG),
template_poses_(template_poses),
align_option_(align_option),
align_chunk_(RG_),
max_registry_shift_input_(max_registry_shift)
{
	bool alignment_from_template = option[challenge::aligned]();
	
	// set up secstruct chunks
	template_ss_chunks_.clear();
	template_ss_chunks_.resize(template_poses_.size());
	for (core::Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		// find ss chunks in template
		template_ss_chunks_[i_template] = extract_secondary_structure_chunks(*template_poses_[i_template]);
	}
	
	Size count = 0;
	for (core::Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		if (template_ss_chunks_[i_template].size() != 0) ++count;
	}
	if (count == 0) {
		utility_exit_with_message("Template structures need at least one secondary structure for this protocol");
	}
	
	sequence_alignments_.clear();
	for (core::Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		std::map <core::Size, core::Size> sequence_alignment;
		if (alignment_from_template) {
			get_alignment_from_template(template_poses_[i_template], sequence_alignment);
		}
		else {
			std::map <core::Size, core::Size> chunk_mapping;
			if (option[challenge::chunk_mapping].user()) {
				for (Size i=1;i<=option[challenge::chunk_mapping]().size();++i) {
					chunk_mapping[i] = option[challenge::chunk_mapping]()[i];
				}
			}
			
			get_alignment_from_chunk_mapping(chunk_mapping, template_ss_chunks_[i_template], ss_chunks_pose, sequence_alignment);
		}
		sequence_alignments_.push_back(sequence_alignment);
	}
}

void
FoldTreeHybridize::get_alignment_from_template(core::pose::PoseCOP const template_pose, std::map <core::Size, core::Size> & seqpos_alignment) {
	// specific to this case, alignment comes from residue number
	for (core::Size ires=1; ires<=template_pose->total_residue(); ++ires) {
		seqpos_alignment[template_pose->pdb_info()->number(ires)] = ires;
	}
}

	
void
FoldTreeHybridize::get_alignment_from_chunk_mapping(std::map <core::Size, core::Size> const & chunk_mapping,
								 Loops const template_ss_chunks,
								 Loops const target_ss_chunks,
								 std::map <core::Size, core::Size> & sequence_alignment)
{
	max_registry_shift_.resize(target_ss_chunks.size());
	for (Size i_chunk_pose = 1; i_chunk_pose <= target_ss_chunks.size(); ++i_chunk_pose) {
		if (chunk_mapping.find(i_chunk_pose) == chunk_mapping.end()) continue;
		core::Size j_chunk_template = chunk_mapping.find(i_chunk_pose)->second;
		
		Size respos_mid_pose = (target_ss_chunks[i_chunk_pose].start() + target_ss_chunks[i_chunk_pose].stop()) / 2;
		Size respos_mid_template = (template_ss_chunks[j_chunk_template].start() + template_ss_chunks[j_chunk_template].stop()) / 2;
		int offset = respos_mid_template - respos_mid_pose;
		
		using namespace ObjexxFCL::fmt;
		if (target_ss_chunks[i_chunk_pose].length() <= template_ss_chunks[j_chunk_template].length()) {
			max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + template_ss_chunks[j_chunk_template].length() - target_ss_chunks[i_chunk_pose].length();
			for (Size ires=target_ss_chunks[i_chunk_pose].start(); ires<=target_ss_chunks[i_chunk_pose].stop(); ++ires) {
				sequence_alignment[ires] = ires+offset;
				//std::cout << I(4, ires) << I(4, ires+offset) << std::endl;
			}
		}
		else {
			max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + target_ss_chunks[i_chunk_pose].length() - template_ss_chunks[j_chunk_template].length();
			for (Size ires_templ=template_ss_chunks[j_chunk_template].start(); ires_templ<=template_ss_chunks[j_chunk_template].stop(); ++ires_templ) {
				sequence_alignment[ires_templ-offset] = ires_templ;
				//std::cout << I(4, ires_templ) << I(4, ires_templ-offset) << std::endl;
			}
		}
	}
}
	
void FoldTreeHybridize::pick_random_template() {
	assert(template_poses_.size() != 0);
	
	template_number_ = 0;
	while (!template_number_) {
		template_number_ = RG_.random_range(1, template_poses_.size());
		if (template_ss_chunks_[template_number_].size() == 0) template_number_ = 0;
	}
}

void FoldTreeHybridize::pick_random_chunk(core::pose::Pose & pose) {
	jump_number_ = RG_.random_range(1, pose.num_jump());
}

Size FoldTreeHybridize::trial_counter(Size ires) {
	return align_chunk_.trial_counter(ires);	
}
	
void
FoldTreeHybridize::apply(core::pose::Pose & pose) {
	max_registry_shift_.resize(pose.num_jump(), max_registry_shift_input_);

	// pick a random template
	pick_random_template();
	align_chunk_.set_template(template_poses_[template_number_], sequence_alignments_[template_number_]);
	
	// random chunk or loop all chunks
	if (align_option_ == random_chunk) {
		// pick a random jump
		pick_random_chunk(pose);
		align_chunk_.set_aligned_chunk(pose, jump_number_);
		align_chunk_.set_reset_torsion_unaligned(false);
		
		// apply alignment
		int registry_shift = RG_.random_range(-max_registry_shift_[jump_number_], max_registry_shift_[jump_number_]);
		align_chunk_.set_registry_shift(registry_shift);
		align_chunk_.apply(pose);
	}
	else {
		// loop over all jumps
		for (core::Size jump_number=1; jump_number<=pose.num_jump(); ++jump_number) {
			align_chunk_.set_aligned_chunk(pose, jump_number);
			
			// apply alignment
			int registry_shift = RG_.random_range(-max_registry_shift_[jump_number], max_registry_shift_[jump_number]);
			align_chunk_.set_registry_shift(registry_shift);
			align_chunk_.apply(pose);
		}
	}
}

std::string
	FoldTreeHybridize::get_name() const {
	return "FoldTreeHybridize";
}
	
private:
	FoldTreeHybridize align_chunk_;
	numeric::random::RandomGenerator & RG_;
	AlignOption align_option_;
	
	utility::vector1 < core::pose::PoseCOP > template_poses_;
	utility::vector1 < Loops > template_ss_chunks_;
	utility::vector1 < std::map <core::Size, core::Size> > sequence_alignments_;
	Size max_registry_shift_input_;
	utility::vector1 < Size > max_registry_shift_;
	
	Size template_number_; // the jump to be realigned
	Size jump_number_; // the jump to be realigned
}; //class FoldTreeHybridize

} // hybridize 
} // comparative_modeling 
} // protocols

#endif
