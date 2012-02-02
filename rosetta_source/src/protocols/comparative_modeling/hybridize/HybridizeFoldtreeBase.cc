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
/// @file protocols/comparative_modeling/hybridize/HybridizeFoldtreeBase.cc
/// @author Yifan Song

// Unit headers
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
#include <core/util/kinematics_util.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
namespace protocols {
namespace comparative_modeling {
namespace hybridize {

HybridizeFoldtreeBase::HybridizeFoldtreeBase() : virtual_res_(-1) {}

void HybridizeFoldtreeBase::save_foldtree(core::pose::Pose const & pose) {
	saved_ft_ = pose.conformation().fold_tree();
	saved_n_residue_ = pose.total_residue();
}
	
void HybridizeFoldtreeBase::restore_foldtree(core::pose::Pose & pose) {
	if (pose.total_residue() > saved_n_residue_) {
		pose.conformation().delete_residue_range_slow(saved_n_residue_+1, pose.total_residue());
	}
	pose.conformation().fold_tree( saved_ft_ );
}


void HybridizeFoldtreeBase::set_chunks(
	const protocols::loops::Loops & chunks,
	const utility::vector1 < core::Size > & anchor_positions )
{
	chunks_last_ = chunks_;
	anchor_positions_last_ = anchor_positions_;
	
	chunks_ = chunks;
	anchor_positions_ = anchor_positions;
}

void HybridizeFoldtreeBase::update(core::pose::Pose & pose) {
	using core::Size;
	using core::Real;
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	using utility::vector1;
	
	assert(chunks_.num_loop());
	
	// Define jumps, cuts
	vector1<int> cuts;
	vector1<std::pair<int, int> > jumps;
	for (Size i = 1; i <= chunks_.num_loop(); ++i) {
		const Loop& chunk = chunks_[i];
		const Size cut_point  = chunk.stop();
		const Size jump_point = anchor_positions_[i];
		
		cuts.push_back(cut_point);
		jumps.push_back(std::make_pair(virtual_res_, jump_point));
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
	
	// Construct the star fold tree from the set of jumps and cuts above.
	// Reorder the resulting fold tree so that <virtual_res> is the root.
	core::kinematics::FoldTree tree(pose.fold_tree());
	bool status = tree.tree_from_jumps_and_cuts(virtual_res_,   // nres_in
												jumps.size(),   // num_jump_in
												ft_jumps,	   // jump_point
												ft_cuts,		// cuts
												virtual_res_);  // root
	if (!status) {
		utility_exit_with_message("HybridizeFoldtreeBase: failed to build fold tree from cuts and jumps");
	}
	
	// Update the pose's fold tree
	pose.fold_tree(tree);
}

	
/// Note: assumes <chunks> are sorted in increasing order of start position
void HybridizeFoldtreeBase::initialize(core::pose::Pose & pose) {
  using core::Size;
  using core::Real;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using utility::vector1;

  assert(chunks_.num_loop());

  // Number of residues before addition of virtual residue
  num_nonvirt_residues_ = pose.total_residue();

  // Add a virtual residue at the center of mass (updates the fold tree)
  numeric::xyzVector<Real> center;
  core::pose::addVirtualResAsRoot(center, pose);

  // Initialize member variable <virtual_res_> with the index of the newly added
  // virtual residue. Subsequently, <virtual_res_> can serve as a proxy for
  // <num_residues>
  virtual_res_ = pose.total_residue();

	update(pose);
	core::util::add_cutpoint_variants(&pose);
}


}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols
