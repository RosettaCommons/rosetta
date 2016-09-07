// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CommonFragmentComparators.hh
/// @brief provides a few ways to compare fragment candidates
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CommonFragmentComparators_hh
#define INCLUDED_protocols_frag_picker_CommonFragmentComparators_hh

// package headers
#include <protocols/frag_picker/CommonFragmentComparators.fwd.hh>
#include <protocols/frag_picker/FragmentComparatorBase.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {


/// @brief comparator based on a sequence position in a query
class CompareQueryPosition : public FragmentComparatorBase {

public:

	/// @brief compares two fragment candidates
	/// @details returns true if second_candidate starts on the later residue than first_candidate
	/// Using this comparator to sort fragment candidates orders them ascending by their sequence position
	bool operator()(
		ScoredCandidate first_candidate,
		ScoredCandidate second_candidate) override {
		return (first_candidate.first->get_first_index_in_query()
			< second_candidate.first->get_first_index_in_query());
	}

};

/// @brief comparator based on the total score of fragments
/// @details returns true if first pair has lower total score than the second one
class CompareTotalScore : public FragmentComparatorBase {

public:
	/// @brief Sets up the comparator for a given FragmentScoreManager.
	/// @details Only the scoring manager knows how to calculate the total score
	/// from a vector of small scores (those gathered in a FragmentScoreMap object)
	CompareTotalScore(scores::FragmentScoreManagerOP scoring) {
		scoring_ = scoring;
	}

	/// @brief compares two fragment candidates
	/// @details returns true if second pair has greater total score than the first one.
	/// Using this comparator to sort fragment candidates order them descending according
	/// to their total score
	bool operator()(
		ScoredCandidate first_candidate,
		ScoredCandidate second_candidate) override {
		return (scoring_->total_score(first_candidate.second)
			< scoring_->total_score(second_candidate.second));
	}

	~CompareTotalScore() override = default;
private:
	scores::FragmentScoreManagerOP scoring_;
};

/// @brief comparator based on the linear combination of some score components
/// @details returns true if first pair has lower weighted score than the second one
class CompareByScoreCombination : public FragmentComparatorBase {

public:
	/// @brief Sets up the comparator for a given FragmentScoreManager.
	/// @details Only the scoring manager knows how to calculate the total score
	/// from a vector of small scores (those gathered in a FragmentScoreMap object)
	CompareByScoreCombination(utility::vector1<Size> which_components,utility::vector1<Real> weights) {

		debug_assert ( which_components.size() == weights.size() );
		for ( Size i=1; i<=which_components.size(); i++ ) {
			components_.push_back( which_components[i] );
			weights_.push_back( weights[i] );
		}
	}

	/// @brief compares two fragment candidates
	/// @details returns true if second pair has greater total score than the first one.
	/// Using this comparator to sort fragment candidates order them descending according
	/// to their total score
	bool operator()(
		ScoredCandidate first_candidate,
		ScoredCandidate second_candidate) override {

		Real t1(0);
		Real t2(0);
		for ( Size i=1; i<=components_.size(); i++ ) {
			t1 += first_candidate.second->at( components_[i] ) * weights_[i];
			t2 += second_candidate.second->at( components_[i] ) * weights_[i];
		}
		return (t1 < t2);
	}

	~CompareByScoreCombination() override = default;
private:
	utility::vector1<Size> components_;
	utility::vector1<Real> weights_;
};


/// @brief comparator based on one of the score components calculated for  fragments
/// @details returns true if first pair has lower score component than the second one
class CompareScoreComponent : public FragmentComparatorBase {

public:
	/// @brief Sets up the comparator based on a given score component
	CompareScoreComponent(Size component_id) {
		component_id_ = component_id;
	}

	/// @brief compares two fragment candidates
	/// @details returns true if second pair has greater total score than the first one.
	/// Using this comparator to sort fragment candidates order them ascending according
	/// to their total score
	bool operator()(
		ScoredCandidate first_candidate,
		ScoredCandidate second_candidate) override {
		return (first_candidate.second->get_score_components()[component_id_]
			< second_candidate.second->get_score_components()[component_id_]);
	}

private:
	Size component_id_;
};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_CommonFragmentComparators_HH */
