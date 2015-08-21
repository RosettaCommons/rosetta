// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.hh
/// @brief  scores a fragment by Blosum62 profile score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreBlosum62_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreBlosum62_hh

// type headers
#include <core/types.hh>


// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  a fragment candidate
class ProfileScoreBlosum62: public CachingScoringMethod {
public:

	ProfileScoreBlosum62(Size priority, Real lowest_acceptable_value, bool use_lowest,
		sequence::SequenceProfileOP query_profile, Size longest_vall_chunk);

	~ProfileScoreBlosum62();

	void do_caching(VallChunkOP);
	void clean_up() {
	}
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	//bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);
protected:
	Matrix scores_;

private:
	sequence::SequenceProfileOP query_profile_;
	utility::vector1< utility::vector1< Real > > blosum_matrix_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScoreBlosum62: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreBlosum62() :
		MakeFragmentScoringMethod("ProfileScoreBlosum62") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreBlosum62_HH */
