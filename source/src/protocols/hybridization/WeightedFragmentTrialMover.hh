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

#ifndef INCLUDED_protocols_hybridization_WeightedFragmentTrialMover_hh
#define INCLUDED_protocols_hybridization_WeightedFragmentTrialMover_hh

#include <protocols/hybridization/WeightedFragmentTrialMover.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <core/fragment/FragSet.hh>

#include <numeric/random/WeightedSampler.hh>

namespace protocols {
namespace hybridization {

class WeightedFragmentTrialMover : public protocols::moves::Mover
{
public:
	WeightedFragmentTrialMover(
		utility::vector1< core::fragment::FragSetOP > const & frag_libs,
		utility::vector1< core::Real > const & residue_weights,
		utility::vector1< core::Size > const & anchor_residues=utility::vector1< core::Size >(0),
		core::Size const nr_frags = 0
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
};

} // hybridization
} // protocols

#endif
