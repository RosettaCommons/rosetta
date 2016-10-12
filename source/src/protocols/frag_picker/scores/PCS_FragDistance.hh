// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/MidPhiOut.hh
/// @brief  Scores fragments by disulfide-linke Calpha distances
/// @author Robert Vernon

#ifndef INCLUDED_protocols_frag_picker_scores_PCS_FragDistance_HH
#define INCLUDED_protocols_frag_picker_scores_PCS_FragDistance_HH

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>

#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini
#include <core/types.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  scores a fragment by the root mean square deviation of Phi and Psi angles.
class PCS_FragDistance: public CachingScoringMethod {
public:

	/// @brief  creates a Calpha distance based scoring function.
	/// @details Scores distance matches for short range disulfides (ie: within fragment)
	PCS_FragDistance(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		ObjexxFCL::FArray2D_double target_ca_dev,
		ObjexxFCL::FArray2D_double target_ca_dist, core::Size largest_fragment, core::Size max_length);

	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	core::Size n_res_;
	std::string cached_scores_id_;

	ObjexxFCL::FArray2D_double chunk_ca_distances_;
	ObjexxFCL::FArray2D_double target_ca_dev_;
	ObjexxFCL::FArray2D_double target_ca_dist_;

	core::Size largest_fragment_;
	core::Size max_length_;
};

/// @brief  Matker class that produces a new PCS_FragDistance object
class MakePCS_FragDistance: public MakeFragmentScoringMethod {
public:

	MakePCS_FragDistance() :
		MakeFragmentScoringMethod("PCS_FragDistance") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_PCS_FragDistance_HH
