// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentScoreManager.hh
/// @brief  Score manager
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentScoreManager_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentScoreManager_hh

// package headers
#include <protocols/frag_picker/scores/FragmentScoreManager.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++
#include <map>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief holds particular score components, weights and calculates the total score for a fragment candidate
/// @details a fragment picker object needs exactly one fragment manager to pick fragments. Adding new scoring methods
/// is supposed to be done FragmentPicker, which calls proper method from this class.
class FragmentScoreManager: public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FragmentScoreManager();

	/// @brief precision used to display the total score for each fragment
	static const Size TOTAL_PRECISION = 3;
	/// @brief default width used to print a score value, equal to 6
	Size default_width_;
	/// @brief default precision used to print a score value, equal to 1
	Size default_precision_;

	/// @brief creates an empty manager
	FragmentScoreManager();

	/// @brief creates an empty score map
	FragmentScoreMapOP create_empty_map();

	/// @brief registers a new scoring method in this manager
	void add_scoring_method(
		FragmentScoringMethodOP, Real);

	/// @brief says how many scoring methods have already been registered
	inline Size count_components() {
		return scores_.size();
	}

	/// @brief Returns a desired scoring method
	/// @details Allowed index values are [1,count_components()]
	inline FragmentScoringMethodOP get_component(
		Size index) {
		return scores_[index];
	}

	/// @brief returns a vector of weights that are used to compute the total score
	utility::vector1<Real> get_weights() {
		return score_weights_;
	}

	/// @brief prints a nice table showing the registered scores
	/// @details the table shows also the order in which the scores are evaluated
	void show_scoring_methods(std::ostream& out);

	/// @brief calculates the total score
	Real total_score(FragmentScoreMapOP);

	/// @brief calculates all the small scores for a given fragment
	/// @brief results are properly stored inside a FragmentScoreMap object
	virtual bool score_fragment_from_cache(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief calculates all the small scores for a given fragment
	/// @brief results are properly stored inside a FragmentScoreMap object
	virtual bool score_fragment(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief those score metohods that have weight = 0.0 will be computed after the fragments are picked
	/// @details if weight for a score is 0.0 than it does not affect the total score and thus has no
	///  effect on fragments sorting, quota, etc. Such scores may be computed after fragment picking is finished
	// By default it is set to false.
	void use_late_scoring_for_zeros(const bool if_true) {
		zeros_score_later_ = if_true;
	}

	/// @brief says if late scoring is used or not
	/// @details late scoring means that some scores (those with weight=0, such as crmsd) are evaluated
	/// only for the selected fragments rather than for all the candidates
	inline bool if_late_scoring_for_zeros() {
		return zeros_score_later_;
	}

	/// @brief calculates all these small scores for a given fragment whose weight is 0.0
	/// @details When use_late_scoring_for_zeros() was used to set the flag to true,
	/// all fragment scoring methods neglects zero-weighted scores.
	/// These will be evaluated by this function, that may be called after fragments are picked.
	/// This way some time consuming computations (e.g. crmsd for fragments) may be computed
	/// only for the selected fragments rather than for all of them
	bool score_zero_scores(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief registers a maker object that will be used to create a scoring method object
	void register_score_maker(MakeFragmentScoringMethodOP);

	/// @brief creates a new scoring method object
	void create_scoring_method(std::string const &, Size, Real, Real, bool,
		FragmentPickerOP, std::string);

	/// @brief reads a config file and creates scoring methods
	void create_scores(std::string const &, FragmentPickerOP);

	/// @brief calls do_caching() for each FragmentScoringMethod object, if it is possible
	/// @details FragmentPicker calls this method when a new chunk is to be processed
	void do_caching(VallChunkOP chunk);

	/// @brief calls clean_up() for each FragmentScoringMethod object, if it is possible
	/// @details FragmentPicker calls this method when a given chunk has been processed
	void clean_up();

	/// @brief sets up a new number of characters spend to print fragment score value
	void set_width(FragmentScoringMethodOP which_score, Size width) {
		width_[which_score] = width;
	}

	/// @brief sets up a new precision used to print fragment score value
	void set_precision(FragmentScoringMethodOP which_score, Size precision) {
		precision_[which_score] = precision;
	}

	/// @brief prints a flat table with all scores for all the fragments in a given vector
	virtual void describe_fragments(utility::vector1<std::pair<FragmentCandidateOP,
		scores::FragmentScoreMapOP> > const&, std::ostream&);

protected:
	utility::vector1<FragmentScoringMethodOP> scores_;
	std::map<FragmentScoringMethodOP, Size> width_;
	std::map<FragmentScoringMethodOP, Size> precision_;
private:
	utility::vector1<Real> score_weights_;
	std::map<std::string, MakeFragmentScoringMethodOP> registered_makers_;
	bool zeros_score_later_;
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_FragmentScoreManager_HH */

