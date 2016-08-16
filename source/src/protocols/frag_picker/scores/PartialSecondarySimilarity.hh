// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/PartialSecondarySimilarity.cc
/// @brief  scores a fragment by secondary structure similarity
///         but throws out the worst 20% of the residue matches
/// @author Robert Vernon (rvernon@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_scores_PartialSecondarySimilarity_hh
#define INCLUDED_protocols_frag_picker_scores_PartialSecondarySimilarity_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>

#include <core/fragment/SecondaryStructure.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  scores a fragment by secondary structure similarity
/// @detail The score for each position is P(H), P(L) or P(E) if
/// a vall residue is within Helix, Loop or Extended secondary stucture element, respectively.
/// P(H), P(L) and P(E) denotes the probability that a given residue in a query
/// is within Helix, Loop or Extended secondary stucture element.
/// The total score of a fragment is a simple sum of all positions; for N-mer fragment is a sum of N terms\n
/// If P(H), P(L) and P(E) probabilities takes only 1.0 and 0.0 values, result of this scoring function
/// should be the same as SecondaryIdentity, although the later one is faster.
class PartialSecondarySimilarity: public CachingScoringMethod {
public:

	PartialSecondarySimilarity(Size priority, Real lowest_acceptable_value, bool use_lowest,
		core::fragment::SecondaryStructureOP query_prediction, std::string prediction_name,
		Size sequence_length, Size longest_vall_chunk);

	~PartialSecondarySimilarity() {}

	void do_caching(VallChunkOP);

	bool cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map);
	void clean_up() {}

	/// @brief Computes the score
	virtual bool score(FragmentCandidateOP, FragmentScoreMapOP);


	inline std::string& get_prediction_name() { return prediction_name_; };
protected:
	Matrix scores_;

private:
	std::string prediction_name_;

	utility::vector1< Real > norm_query_H_;
	utility::vector1< Real > norm_query_E_;
	utility::vector1< Real > norm_query_L_;

	utility::vector1< utility::vector1< Real > > raw_probs_;
	Size query_len_;

	std::string cached_scores_id_;
};

/// @brief  Maker class that produces a new PartialSecondarySimilarity object
class MakePartialSecondarySimilarity: public MakeFragmentScoringMethod {
public:

	MakePartialSecondarySimilarity() :
		MakeFragmentScoringMethod("PartialSecondarySimilarity") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value, bool use_lowest,
		FragmentPickerOP picker, std::string prediction_id) {

		Size sequence_length = picker->get_query_seq()->length();
		Size vall_max_len = picker->get_vall()->get_largest_chunk_size();

		core::fragment::SecondaryStructureOP query_prediction( picker->get_query_ss(prediction_id) );
		if ( ! query_prediction ) {
			utility_exit_with_message("Unable to find secondary structure prediction for "+prediction_id);
		}
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PartialSecondarySimilarity(priority,
			lowest_acceptable_value, use_lowest,
			query_prediction,prediction_id,
			sequence_length,vall_max_len) );
	}
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_PartialSecondarySimilarity_HH */
