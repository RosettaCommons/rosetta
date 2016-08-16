// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScoreStructL1.hh
/// @brief  scores a fragment by L1 structure based profile score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)
/// @author David E Kim (dekim@uw.edu)

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreStructL1_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreStructL1_hh

// type headers
#include <core/types.hh>


// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/fragment_scoring_utilities.hh>

#include <core/sequence/SequenceProfile.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  a fragment candidate
class ProfileScoreStructL1: public CachingScoringMethod {
public:

	ProfileScoreStructL1(Size, Real , bool,
		sequence::SequenceProfileOP, utility::vector1<Size> &,Size);
	~ProfileScoreStructL1();

	void do_caching_simple(VallChunkOP);
	void do_caching(VallChunkOP);
	void clean_up() {
	}
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	//bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);
protected:
	Matrix scores_;
	utility::vector1< Matrix > cache_;

private:
	sequence::SequenceProfileOP query_profile_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScoreStructL1: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreStructL1() :
		MakeFragmentScoringMethod("ProfileScoreStructL1") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreStructL1_HH */
