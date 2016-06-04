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
/// @author David Kim, Yifan Song

#ifndef INCLUDED_protocols_hybridization_WeightedFragmentSmoothTrialMover_hh
#define INCLUDED_protocols_hybridization_WeightedFragmentSmoothTrialMover_hh

#include <protocols/hybridization/WeightedFragmentSmoothTrialMover.fwd.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/types.hh>
#include <basic/prof.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <cstdlib>
#include <string>

namespace protocols {
namespace hybridization {

class WeightedFragmentSmoothTrialMover : public protocols::moves::Mover
{
public:
	WeightedFragmentSmoothTrialMover(
		utility::vector1< core::fragment::FragSetOP > const & frag_libs,
		utility::vector1< core::Real > const & residue_weights,
		utility::vector1< core::Size > const & anchor_residues,
		core::Size const nr_frags,
		simple_moves::FragmentCostOP cost
	);
	void update_sampler_weights( utility::vector1< core::Real > const residue_weights );

	void apply(core::pose::Pose & pose);
	std::string get_name() const;
	core::Size get_total_frames() { return total_frames_; }
	core::Size get_nr_frags() { return nr_frags_; }
private:
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
	utility::vector1< numeric::random::WeightedSampler > weighted_sampler_;
	utility::vector1< Size > anchor_reses_; // list of residue indices
	core::Size nr_frags_; // number of fragments to use per frame (1 to this number), if 0, use all
	core::Size total_frames_;
	simple_moves::FragmentCostOP cost_;

	// choose randomly fragments that are below cutoff_
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real cutoff_;
};

} // hybridization
} // protocols

#endif
