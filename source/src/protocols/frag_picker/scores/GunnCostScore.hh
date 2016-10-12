// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/GunnCostScore.hh
/// @brief  Object that scores a fragment by Gunn cost function (how it fits to the local structure)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_GunnCostScore_hh
#define INCLUDED_protocols_frag_picker_scores_GunnCostScore_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/GunnCost.hh>


#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<core::Real> > Matrix;

/// @brief  scores a fragment by its crmsd to the given reference structure
class GunnCostScore: public CachingScoringMethod {
public:

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to a given pose, which should have the same number of residues a the query sequence
	GunnCostScore(core::Size, core::Real, bool, core::pose::PoseOP,utility::vector1<core::Size>,core::Size);

	~GunnCostScore();

	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

	void computeGunnTuples(core::pose::Pose &,core::Size,utility::vector1<GunnTuple> &);

private:
	core::Size n_atoms_;
	core::pose::PoseOP reference_pose_;
	utility::vector1<core::Size> frag_sizes_;
	core::Size max_chunk_size_;
	utility::vector1<utility::vector1<GunnTuple> > ref_gunn_data_;
	utility::vector1<utility::vector1<GunnTuple> > chunk_gunn_data_;
	GunnCost gunn_cost_;
};

/// @brief  Maker class that produces a new GunnCostScore object
class MakeGunnCostScore: public MakeFragmentScoringMethod {
public:

	MakeGunnCostScore() :
		MakeFragmentScoringMethod("GunnCostScore") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_GunnCostScore_HH
