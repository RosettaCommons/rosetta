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
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/util/kinematics_util.hh>


//Auto Headers
#include <core/conformation/Conformation.hh>
namespace protocols {
namespace comparative_modeling {
namespace hybridize {

HybridizeFoldtreeDynamic::HybridizeFoldtreeDynamic() {}
	
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
	
void HybridizeFoldtreeDynamic::make_complete_chunks(
	utility::vector1 < core::Size > cut_positions,
	core::Size n_residues
)
{
	assert(cut_positions.size() == core_chunks_.size() - 1);

	chunks_ = core_chunks_;
	for (core::Size i=1; i<=chunks_.num_loop(); ++i) {
			if (i==1) {
				chunks_[i].set_start(1);
			} else {
				core::Size new_start = cut_positions[i-1] + 1;
				chunks_[i].set_start( new_start );
			}

			if (i==chunks_.num_loop()) {
				chunks_[i].set_stop(n_residues);
			} else {
				core::Size new_stop = cut_positions[i];
				chunks_[i].set_stop( new_stop );
			}
	}
}

void HybridizeFoldtreeDynamic::choose_anchors() {
	anchor_positions_.clear();
	for (core::Size i=1; i<=core_chunks_.num_loop(); ++i) {
		anchor_positions_.push_back( choose_anchor_position(core_chunks_[i]) );
	}		
}

/// from cmiles:
///   mu : midpoint of the chunk
///   sigma : linear function of chunk length
core::Size HybridizeFoldtreeDynamic::choose_anchor_position(const protocols::loops::Loop & chunk) const {
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
  using core::Size;
  using core::Real;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using utility::vector1;

	num_nonvirt_residues_ = pose.total_residue();

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		num_nonvirt_residues_ = symm_info->num_independent_residues();
	}
	while (!pose.residue(num_nonvirt_residues_).is_protein()) num_nonvirt_residues_--;

	set_core_chunks(core_chunks);

	// save previous
	chunks_last_ = chunks_;
	anchor_positions_last_ = anchor_positions_;

	// set new
	choose_anchors();
	utility::vector1 < core::Size > cut_positions = decide_cuts(num_nonvirt_residues_);
	make_complete_chunks(cut_positions, num_nonvirt_residues_);

  assert(chunks_.num_loop());

  // Add a virtual residue at the center of mass (updates the fold tree)
	if (!symm_info)  // pose is asymm
	  core::pose::addVirtualResAsRoot(pose);
  virtual_res_ = pose.total_residue();

	// do the actual foldtree updates
	update(pose);

	protocols::loops::add_cutpoint_variants( pose );
}


void HybridizeFoldtreeDynamic::update(core::pose::Pose & pose) {
	using core::Size;
	using core::Real;
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	using utility::vector1;
	
	assert(chunks_.num_loop());

	bool use_symm = core::pose::symmetry::is_symmetric(pose);

	// Define jumps, cuts
	vector1<int> cuts;
	vector1<std::pair<int, int> > jumps;
	for (Size i = 1; i <= chunks_.num_loop(); ++i) {
		const Loop& chunk = chunks_[i];
		const Size cut_point  = chunk.stop();
		const Size jump_point = anchor_positions_[i];

		Size jump_root = num_nonvirt_residues_+1;
		if (use_symm && i>1) jump_root = anchor_positions_[1];
		
		cuts.push_back(cut_point);
		jumps.push_back(std::make_pair(jump_root, jump_point));
	}
	
	// Remember to include the original cutpoint at the end of the chain
	// (before the virtual residue)
	cuts.push_back(num_nonvirt_residues_);
	
	ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
	for (Size i = 1; i <= jumps.size(); ++i) {
		ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
		ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
	}
	
	ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
	for (Size i = 1; i <= cuts.size(); ++i) {
		ft_cuts(i) = cuts[i];
	}
	
	// "symmetry-safe" version
	core::kinematics::FoldTree tree = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	bool status = tree.tree_from_jumps_and_cuts(num_nonvirt_residues_+1,   // nres_in
												jumps.size(),   // num_jump_in
												ft_jumps,	   // jump_point
												ft_cuts,		// cuts
												num_nonvirt_residues_+1);  // root
	if (!status) {
		utility_exit_with_message("HybridizeFoldtreeBase: failed to build fold tree from cuts and jumps");
	}

	std::cerr << tree;
	
	// Update the pose's fold tree
	core::pose::symmetry::set_asymm_unit_fold_tree( pose , tree );
}

	
void HybridizeFoldtreeDynamic::set_core_chunks(const protocols::loops::Loops & chunks) {
	core_chunks_last_ = core_chunks_;
	core_chunks_ = chunks;
}
	
}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols
