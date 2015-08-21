// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/StarTreeBuilder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.hh>

// C/C++ headers
#include <string>
#include <utility>

// External headers
#include <boost/unordered/unordered_map.hpp>

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
namespace nonlocal {

// Static member initialization
const std::string StarTreeBuilder::PREFIX_INITIAL = "initial";
const std::string StarTreeBuilder::PREFIX_FINAL = "final";

StarTreeBuilder::StarTreeBuilder() : virtual_res_(-1) {}

/// Note: assumes <chunks> are sorted in increasing order of start position
void StarTreeBuilder::set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) {
	using core::Size;
	using core::Real;
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	using utility::vector1;

	assert(pose);
	assert(chunks.num_loop());

	// Number of residues before addition of virtual residue
	Size num_residues = pose->total_residue();

	// Add a virtual residue at the center of mass (updates the fold tree)
	numeric::xyzVector<Real> center;
	chunks.center_of_mass(*pose, &center);
	core::pose::addVirtualResAsRoot(center, *pose);

	// Initialize member variable <virtual_res_> with the index of the newly added
	// virtual residue. Subsequently, <virtual_res_> can serve as a proxy for
	// <num_residues>
	virtual_res_ = pose->total_residue();

	// Define jumps, cuts
	vector1<int> cuts;
	vector1<std::pair<int, int> > jumps;
	for ( Loops::const_iterator i = chunks.begin(); i != chunks.end(); ++i ) {
		const Loop& chunk = *i;
		const Size cut_point  = chunk.stop();
		const Size jump_point = choose_anchor_position(chunk);

		cuts.push_back(cut_point);
		jumps.push_back(std::make_pair(virtual_res_, jump_point));
	}

	// Remember to include the original cutpoint at the end of the chain
	// (before the virtual residue)
	cuts.push_back(num_residues);

	ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
	for ( Size i = 1; i <= jumps.size(); ++i ) {
		ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
		ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
	}

	ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		ft_cuts(i) = cuts[i];
	}

	// Construct the star fold tree from the set of jumps and cuts above.
	// Reorder the resulting fold tree so that <virtual_res> is the root.
	core::kinematics::FoldTree tree(pose->fold_tree());
	bool status = tree.tree_from_jumps_and_cuts(virtual_res_,   // nres_in
		jumps.size(),   // num_jump_in
		ft_jumps,       // jump_point
		ft_cuts,        // cuts
		virtual_res_);  // root
	if ( !status ) {
		utility_exit_with_message("StarTreeBuilder: failed to build fold tree from cuts and jumps");
	}

	// Update the pose's fold tree
	pose->fold_tree(tree);

	// When native is available, compute RMSD over jump residues. Results stored
	// as comments in the pose with format: rmsd_jump_residue_X rmsd
	do_compute_jump_rmsd(pose, PREFIX_INITIAL);
}

void StarTreeBuilder::do_compute_jump_rmsd(core::pose::Pose* model, const std::string& prefix) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using boost::unordered_map;
	using core::Real;
	using core::Size;
	using core::pose::Pose;
	assert(model);

	if ( !option[OptionKeys::in::file::native].user() ) {
		return;
	}

	// Compute RMSD of jump residues (jump point +/- 1 residue)
	unordered_map<Size, Real> rmsds;
	Pose native = *core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
	core::scoring::compute_jump_rmsd(native, *model, &rmsds);

	// Write results to <model> as comments
	for ( unordered_map<Size, Real>::const_iterator i = rmsds.begin(); i != rmsds.end(); ++i ) {
		Size residue = i->first;
		Real rmsd = i->second;
		core::pose::add_comment(*model,
			(boost::format("%1% rmsd_jump_residue_%2%") % prefix % residue).str(),
			(boost::format("%1%") % rmsd).str());
	}
}

void StarTreeBuilder::tear_down(core::pose::Pose* pose) {
	assert(pose);
	if ( virtual_res_ == -1 ) {
		return;
	} else {
		do_compute_jump_rmsd(pose, PREFIX_FINAL);
		pose->conformation().delete_residue_slow(virtual_res_);
		virtual_res_ = -1;
	}
}

/// mu- midpoint of the chunk
/// sigma- linear function of chunk length
core::Size StarTreeBuilder::choose_anchor_position(const protocols::loops::Loop& chunk) const {
	using boost::math::normal;
	using core::Size;
	using numeric::random::DistributionSampler;

	double mu = chunk.start() + chunk.length() / 2.0;
	double sigma = chunk.length() / 5.0;
	normal distribution(mu, sigma);
	DistributionSampler<normal> sampler(distribution);

	// Clamp insertion position to closed interval [start, stop]
	Size position = static_cast<Size>(sampler.sample());
	return numeric::clamp<Size>(position, chunk.start(), chunk.stop());;
}

}  // namespace nonlocal
}  // namespace protocols
