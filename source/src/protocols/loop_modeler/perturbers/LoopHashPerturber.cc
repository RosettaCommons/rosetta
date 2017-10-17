// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// @author Xingjie Pan (xingjiepan@gmail.com)

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/loop_modeler/perturbers/LoopHashPerturber.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loophash/LoopHashLibrary.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.loop_modeler.perturbers.LoopHashPerturber" );

namespace protocols {
namespace loop_modeler {
namespace perturbers {



LoopHashPerturber::LoopHashPerturber(protocols::loophash::LoopHashLibraryOP lh_library){
	lh_library_ = lh_library;
}

void LoopHashPerturber::perturb_subset(
	core::pose::Pose const &pose, kinematic_closure::IndexList const & residues, kinematic_closure::ClosureProblemOP problem) {
	using namespace protocols::loophash;
	using numeric::conversions::DEGREES;

	// Get the backbone segments
	// Because the torsions at the terminal pivots should also be variable, the
	// segement to get should extend beyond the terminal pivots by one residue
	if ( random_mode_ ) {
		get_random_backbone_segments(pose, residues.front() - 1, residues.back() + 1);
	} else {
		get_backbone_segments(pose, residues.front() - 1, residues.back() + 1);
	}

	//TR << "Find " << bb_segs_.size() << " segments" << std::endl;///DEBUG

	if ( bb_segs_.size() == 0 ) {
		return;
	}

	// Pick a random segment
	core::Size seg_id = numeric::random::random_range(1, bb_segs_.size());
	BackboneSegment bb_seg = bb_segs_[seg_id].first;

	// Get the torsions, NOTE: these are std::vectors
	std::vector<core::Real> phis = bb_seg.phi();
	std::vector<core::Real> psis = bb_seg.psi();
	std::vector<core::Real> omegas = bb_seg.omega();

	core::Size start = random_mode_ ? numeric::random::random_range(1, residues.size()) : 1;
	core::Size stop = random_mode_ ? numeric::random::random_range(start, residues.size()) : residues.size();

	// Apply the torsions of the selected segment
	for ( core::Size i=start; i<=stop; ++i ) {

		//TR << "residues.size() = " << residues.size() << ", start = " << start << ", stop = " << stop << ", phi = " << phis[i] << ", psi = " << psis[i] << ", omega = " << omegas[i] << std::endl; ///DEBUG

		problem->perturb_phi(residues[i], phis[i] , DEGREES);
		problem->perturb_psi(residues[i], psis[i] , DEGREES);
		problem->perturb_omega(residues[i], omegas[i] , DEGREES);
	}

	// Perturb the sequence

	if ( perturb_sequence_ ) {
		//std::cout << "sequence: " << bb_segs_[seg_id].second << std::endl;///DEBUG
		utility::vector1 <std::string> sequence = problem->unperturbed_sequence();
		for ( core::Size i=start; i<=stop; ++i ) {
			sequence[i] = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code(bb_segs_[seg_id].second[i]) );
		}
		problem->perturbed_sequence(sequence);
	}
}

void
LoopHashPerturber::get_backbone_segments(core::pose::Pose const& pose,
	core::Size loophash_fragment_start,
	core::Size loophash_fragment_end){
	using namespace protocols::loophash;

	core::Size loop_size = loophash_fragment_end - loophash_fragment_start + 1;

	//Get the loop transform
	numeric::geometry::hashing::Real6 loop_transform;
	if ( !get_rt_over_leap_without_foldtree_bs( pose, loophash_fragment_start, loophash_fragment_end+1, loop_transform ) ) {
		utility_exit_with_message("Unable to find rigid body transform over jump");
	}

	//Update the segments only if the transform is changed considerablely
	core::Real diff2 = 0;
	for ( core::Size i=1; i<=6; ++i ) {
		diff2 += (loop_transform[i] - last_loop_transform_[i]) * (loop_transform[i] - last_loop_transform_[i]);
	}
	if ( diff2 < 1 ) {
		return;
	}
	last_loop_transform_ = loop_transform;

	//Get the hash map
	LoopHashMap & hashmap = lh_library_->gethash(loop_size);

	//Save all backbone segments from this transform
	BackboneSegments bb_segs;

	std::vector<core::Size> leap_index_bucket;

	if ( use_radial_lookup_ ) {
		hashmap.radial_lookup( core::Size(4/*max radius*/), loop_transform, leap_index_bucket);
	} else {
		hashmap.lookup(loop_transform, leap_index_bucket);
	}

	for ( auto const& index : leap_index_bucket ) {
		bb_segs.push_back(extract_fragment(index, loop_size));
	}

	bb_segs_ = bb_segs;
}

void
LoopHashPerturber::get_random_backbone_segments(core::pose::Pose const&,
	core::Size loophash_fragment_start,
	core::Size loophash_fragment_end){

	using namespace protocols::loophash;

	core::Size loop_size = loophash_fragment_end - loophash_fragment_start + 1;
	LoopHashMap & hashmap = lh_library_->gethash(loop_size);

	core::Size frag_index = numeric::random::random_range(0, hashmap.n_loops() - 1);

	bb_segs_.clear();
	bb_segs_.push_back(extract_fragment(frag_index, loop_size));
}

std::pair< protocols::loophash::BackboneSegment, std::string >
LoopHashPerturber::extract_fragment(core::Size frag_index, core::Size loop_size){
	using namespace protocols::loophash;

	LoopHashMap & hashmap = lh_library_->gethash(loop_size);

	BBData bb_data;
	BBExtraData extra_data;
	BackboneSegment backbone_seg;

	// Get the backbone
	LeapIndex cp = hashmap.get_peptide( frag_index );
	lh_library_->backbone_database().get_backbone_segment( cp.index, cp.offset , loop_size , backbone_seg );

	// Get the sequence
	lh_library_->backbone_database().get_protein( cp.index, bb_data );
	lh_library_->backbone_database().get_extra_data(bb_data.extra_key, extra_data);

	runtime_assert(cp.offset % 3 == 0);
	core::Size seq_offset = (cp.offset/3);
	std::string sequence = extra_data.sequence;
	std::string loop_sequence = sequence.substr(seq_offset, loop_size);

	return std::make_pair(backbone_seg, loop_sequence);
}

void LoopHashPerturber::perturb_subset_with_balance(
	core::pose::Pose const & pose, kinematic_closure::IndexList const & residues, kinematic_closure::ClosureProblemOP problem) {

	// TODO: currently the detailed balance is not guaranteed
	perturb_subset(pose, residues, problem);
}

}
}
}
