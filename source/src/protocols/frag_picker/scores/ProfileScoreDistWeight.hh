// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.hh
/// @brief  scores a fragment by weighting L1 profile distances by residue type
/// @author Robert Vernon

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <core/sequence/SequenceProfile.hh>

// utility headers
#include <utility/exit.hh>
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<core::Real> > Matrix;

/// @brief  a fragment candidate
class ProfileScoreDistWeight: public CachingScoringMethod {
public:

	ProfileScoreDistWeight(
		core::Size priority,
		core::Real lowest_acceptable_value,
		bool use_lowest,
		core::sequence::SequenceProfileOP query_profile,
		core::fragment::SecondaryStructureOP query_ss_prediction,
		std::string query_sequence, core::Size longest_vall_chunk
	);

	void do_caching(VallChunkOP) override;
	void clean_up() override {
	}
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP) override;
	//bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map, std::ostream& out);

protected:
	Matrix scores_;

private:
	utility::vector1< utility::vector1< utility::vector1 <core::Real> > > distance_weights_;

	std::string query_sequence_;

	std::map<char,core::Size> ss_type_map_;
	std::map<char,core::Size> aa_order_map_;

	core::sequence::SequenceProfileOP query_profile_;
	core::fragment::SecondaryStructureOP query_ss_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScoreDistWeight: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreDistWeight() :
		MakeFragmentScoringMethod("ProfileScoreDistWeight") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string) override;
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_HH */
