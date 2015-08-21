// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.hh
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScore_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScore_hh

// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScore.fwd.hh>

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/ScoringScheme.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  a fragment candidate
class ProfileScore: public CachingScoringMethod {
public:

	ProfileScore(Size priority, Real lowest_acceptable_value, bool use_lowest,
		sequence::SequenceProfileOP query_profile,
		sequence::ScoringSchemeOP profile_scoring, Size longest_vall_chunk) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "ProfileScore") {
		if ( query_profile->size() != query_profile->length() ) {
			utility_exit_with_message("ProfileScore needs a valid sequence profile.");
		}
		query_profile_ = query_profile;
		profile_scoring_ = profile_scoring;

		for ( Size i = 1; i <= query_profile->length(); ++i ) {
			utility::vector1<Real> row(longest_vall_chunk);
			scores_.push_back(row);
		}
	}

	~ProfileScore();

	void do_caching(VallChunkOP);
	void clean_up() {
	}
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map,
		std::ostream& out);

protected:
	Matrix scores_;

private:
	sequence::SequenceProfileOP query_profile_;
	sequence::ScoringSchemeOP profile_scoring_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScore: public MakeFragmentScoringMethod {
public:

	MakeProfileScore() :
		MakeFragmentScoringMethod("ProfileScore") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScore_HH */
