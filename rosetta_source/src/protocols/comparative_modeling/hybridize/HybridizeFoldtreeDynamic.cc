// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief A mover to dynamically manipulate fold tree during the template hybridization sampling
/// @detailed based on Chris Miles's start tree builder
/// @file protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.cc
/// @author Yifan Song

// Unit headers
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.hh>
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeBase.hh>

// C/C++ headers
#include <string>
#include <utility>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/DistributionSampler.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
namespace protocols {
namespace comparative_modeling {
namespace hybridize {

HybridizeFoldtreeDynamic::HybridizeFoldtreeDynamic()
{}
	
void HybridizeFoldtreeDynamic::choose_anchors() {
	anchor_positions_.clear();
	for (core::Size i=1; i<=core_chunks_.num_loop(); ++i) {
		anchor_positions_.push_back( choose_anchor_position(core_chunks_[i]) );
	}		
}

utility::vector1 < core::Size > HybridizeFoldtreeDynamic::decide_cuts(core::Size n_residues) {
	// complete the chunks to cover the whole protein and customize cutpoints
	// cutpoints in the middle of the loop
	// cut number is the residue before the cut
	std::string cut_point_decision = "middle";

	utility::vector1<bool> cut_options(n_residues, true);

	// not allowing cuts in the core_chunks, currectly is set by the secondary structure of the initial template
	for (core::Size ichunk = 1; ichunk<=core_chunks_.num_loop(); ++ichunk) {
		for (core::Size ires = core_chunks_[ichunk].start(); ires <= core_chunks_[ichunk].stop()-1; ++ires) {
			cut_options[ires] = false;
		}
	}

	utility::vector1 < core::Size > cut_positions;
	for (core::Size i=2; i<=core_chunks_.num_loop(); ++i) {
		core::Size loop_start = core_chunks_[i-1].stop();
		core::Size loop_end = core_chunks_[i].start() - 1;
		core::Size cut;

		if (loop_start >= loop_end) {
			cut = loop_start;
		}
		else {
			if ( cut_point_decision == "middle") {
				cut = (loop_start + loop_end ) /2;
			}
			else if ( cut_point_decision == "combine" ) {
				utility::vector1 < core::Size > cut_residues;
				for (core::Size ires = loop_start; ires <= loop_end; ++ires) {
					if (cut_options[ires]) {
						cut_residues.push_back(ires);
					}
				}
				if (cut_residues.size() > 0) {
					using boost::math::normal;
					using core::Size;
					using numeric::random::DistributionSampler;
					
					double mu = cut_residues.size() / 2.0;
					double sigma = cut_residues.size() / 5.0;
					normal distribution(mu, sigma);
					DistributionSampler<normal> sampler(distribution);
					
					// Clamp insertion position to closed interval [start, stop]
					Size position = static_cast<Size>(sampler.sample());
					position = numeric::clamp<core::Size>(position, 1, cut_residues.size());;
					
					cut = cut_residues[position];
				}
				else {
					cut = (loop_start + loop_end ) /2;
				}
			}
			else {
				utility_exit_with_message("do not know how to make cut points");
			}
		}
		cut_positions.push_back(cut);
	}
	return cut_positions;
}
	
protocols::loops::Loops HybridizeFoldtreeDynamic::make_complete_chunks(
	utility::vector1 < core::Size > cut_positions,
	core::Size n_residues
)
{
	assert(cut_positions.size() == core_chunks_.size() - 1);

	protocols::loops::Loops complete_chunks(core_chunks_); 
	for (core::Size i=1; i<=core_chunks_.num_loop(); ++i) {
			if (i==1) {
				complete_chunks[i].set_start(1);
			}
			else {
				core::Size new_start = cut_positions[i-1] + 1;
				complete_chunks[i].set_start( new_start );
			}
			
			if (i==core_chunks_.num_loop()) {
				complete_chunks[i].set_stop(n_residues);
			}
			else {
				core::Size new_stop = cut_positions[i];
				complete_chunks[i].set_stop( new_stop );
			}
	}
	return complete_chunks;
}
	
/// mu- midpoint of the chunk
/// sigma- linear function of chunk length
/// from cmiles
core::Size HybridizeFoldtreeDynamic::choose_anchor_position(const protocols::loops::Loop & chunk) const
{
	using boost::math::normal;
	using core::Size;
	using numeric::random::DistributionSampler;
	
	double mu = chunk.start() + chunk.length() / 2.0;
	double sigma = chunk.length() / 5.0;
	normal distribution(mu, sigma);
	DistributionSampler<normal> sampler(distribution);
	
	// Clamp insertion position to closed interval [start, stop]
	Size position = static_cast<Size>(sampler.sample());
	return numeric::clamp<core::Size>(position, chunk.start(), chunk.stop());;
}


void HybridizeFoldtreeDynamic::initialize(
	core::pose::Pose & pose,
	protocols::loops::Loops const & core_chunks
) {
	core::Size nres = pose.total_residue();
	while (!pose.residue(nres).is_protein()) nres--;

	set_core_chunks(core_chunks);
	choose_anchors();
	utility::vector1 < core::Size > cut_positions = decide_cuts(pose.total_residue());
	protocols::loops::Loops complete_chunks = make_complete_chunks(cut_positions, nres);
	
	HybridizeFoldtreeBase::set_chunks(complete_chunks, anchor_positions_);
	HybridizeFoldtreeBase::initialize(pose); // initialize foldtree based on complete_chunks
}
	
void HybridizeFoldtreeDynamic::set_core_chunks(const protocols::loops::Loops & chunks) {
	core_chunks_last_ = core_chunks_;
	core_chunks_ = chunks;
}
	
}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols
