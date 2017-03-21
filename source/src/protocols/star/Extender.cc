// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/star/Extender.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/star/Extender.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

namespace protocols {
namespace star {

using core::Real;
using core::Size;
using core::pose::Pose;
using utility::vector1;

static const Real EXT_PHI = -150;
static const Real EXT_PSI = +150;
static const Real EXT_OMG = +180;

void generate_extended_pose(Pose* extended_pose, const std::string& sequence) {
	core::pose::make_pose_from_sequence(*extended_pose, sequence, core::chemical::CENTROID );

	for ( Size i = 1; i <= extended_pose->size(); ++i ) {
		extended_pose->set_phi(i, EXT_PHI);
		extended_pose->set_psi(i, EXT_PSI);
		extended_pose->set_omega(i, EXT_OMG);
	}
}

void copy_residues(const Pose& src, Size start, Size stop, Pose* dst) {
	using core::id::AtomID;
	using core::conformation::Residue;

	for ( Size i = start; i <= stop; ++i ) {
		const Residue& r = src.conformation().residue(i);

		for ( Size j = 1; j <= r.natoms(); ++j ) {
			AtomID id(j, i);
			dst->set_xyz(id, src.xyz(id));
		}
	}
}

Extender::Extender(core::sequence::SequenceAlignmentCOP alignment, int num_residues)
: alignment_(alignment) {
	using core::id::SequenceMapping;
	using protocols::loops::Loops;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	debug_assert(alignment);
	debug_assert(num_residues > 0);

	SequenceMapping mapping = alignment->sequence_mapping(1, 2);
	vector1<int> unaligned_res;

	for ( int i = 1; i < num_residues; ++i ) {
		// unaligned: current residue doesn't map to any residue in the template
		// broke: current residue and next residue are both aligned, but not to consecutive residues in the template
		bool curr_aligned = mapping[i];
		bool next_aligned = mapping[i + 1];
		bool broke = (mapping[i] + 1 != mapping[i + 1]) && next_aligned;

		if ( !curr_aligned ) {
			unaligned_res.push_back(i);
		} else if ( broke ) {
			unaligned_res.push_back(i);
			unaligned_res.push_back(i+1);
		}
	}

	if ( !mapping[num_residues] ) {  // last residue
		unaligned_res.push_back(num_residues);
	}

	std::sort(unaligned_res.begin(), unaligned_res.end());
	vector1<int>::const_iterator i = std::unique(unaligned_res.begin(), unaligned_res.end());
	unaligned_res.resize(i - unaligned_res.begin());

	int min_len = option[OptionKeys::abinitio::star::min_unaligned_len]();
	unaligned_ = protocols::comparative_modeling::pick_loops_unaligned(num_residues, unaligned_res, min_len);
	unaligned_->sequential_order();

	aligned_ = protocols::loops::LoopsOP( new Loops(unaligned_->invert(num_residues)) );
	aligned_->sequential_order();
}

void Extender::extend_unaligned(Pose* pose) {
	using core::kinematics::FoldTree;
	using protocols::loops::Loop;

	if ( aligned_->num_loop() != 2 ) {
		utility_exit_with_message("Unsupported operation");
	}

	// Keep track of the cutpoints we introduce
	cutpoints_.clear();

	const Pose reference = *pose;

	generate_extended_pose(pose, reference.sequence());
	core::scoring::calpha_superimpose_pose(*pose, reference);

	const Loop& f1 = (*aligned_)[1];
	const Loop& f2 = (*aligned_)[2];

	FoldTree l2r, r2l;
	l2r.add_edge(1, pose->size(), core::kinematics::Edge::PEPTIDE);
	r2l.add_edge(pose->size(), 1, core::kinematics::Edge::PEPTIDE);

	pose->fold_tree(r2l);
	copy_residues(reference, f1.start(), f1.stop(), pose);
	core::conformation::idealize_position(f1.start() - 1, pose->conformation());

	pose->fold_tree(l2r);
	core::conformation::idealize_position(f1.stop(), pose->conformation());

	copy_residues(reference, f2.start(), f2.stop(), pose);
	core::conformation::idealize_position(f2.stop(), pose->conformation());

	unsigned cut = choose_cutpoint(f1.stop() + 1, f2.start() - 1);
	cutpoints_.push_back(cut);
}

int Extender::choose_cutpoint(int start, int stop) const {
	using numeric::random::WeightedReservoirSampler;
	debug_assert(start > 0);
	debug_assert(start <= stop);

	if ( !pred_ss_ ) {
		return numeric::random::random_range(start, stop);
	}

	WeightedReservoirSampler<int> sampler(1);

	for ( int i = start; i <= stop; ++i ) {
		double weight = pred_ss_->loop_fraction(i) + 1e-20;
		sampler.consider_sample(i, weight);
	}

	vector1<int> samples;
	sampler.samples(&samples);
	return samples[1];
}

}  // namespace star
}  // namespace protocols
