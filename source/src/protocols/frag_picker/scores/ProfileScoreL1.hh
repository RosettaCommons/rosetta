// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.hh
/// @brief  scores a fragment by L1 profile score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreL1_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreL1_hh

// type headers
#include <core/types.hh>


// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/fragment_scoring_utilities.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

class ProfileScoreL1;
typedef utility::pointer::shared_ptr< ProfileScoreL1 > ProfileScoreL1OP;
typedef utility::pointer::shared_ptr< ProfileScoreL1 const > ProfileScoreL1COP;

/// @brief  a fragment candidate
class ProfileScoreL1: public CachingScoringMethod {
public:

	ProfileScoreL1(Size, Real , bool,
		sequence::SequenceProfileOP, utility::vector1<Size> &,Size);
	~ProfileScoreL1();

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

class MakeProfileScoreL1: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreL1() :
		MakeFragmentScoringMethod("ProfileScoreL1") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreL1_HH */
