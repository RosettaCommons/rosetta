// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScoreSubMatrix.hh
/// @brief  scores a fragment by substitution matrix (e.g. Blosum62) profile score
/// @author Dan Kulp (dwkulp@gmail.com), based on code from Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreSubMatrix_HH
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreSubMatrix_HH

// type headers
#include <core/types.hh>


// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<core::Real> > Matrix;

/// @brief  a fragment candidate
class ProfileScoreSubMatrix: public CachingScoringMethod {
public:

	ProfileScoreSubMatrix(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,              std::string sequence,core::Size longest_vall_chunk,std::string subMatrixFile);
	~ProfileScoreSubMatrix();

	void do_caching(VallChunkOP);
	void clean_up() {
	}
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	// Undefinede, commenting out to fix PyRosetta build  bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);
protected:
	Matrix scores_;

private:
	std::string sequence_;
	std::string subMatrixFile_;

	utility::vector1< utility::vector1< core::Real > > sub_matrix_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScoreSubMatrix: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreSubMatrix() :
		MakeFragmentScoringMethod("ProfileScoreSubMatrix") {
	}

	FragmentScoringMethodOP make(core::Size, core::Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreSubMatrix_HH */
