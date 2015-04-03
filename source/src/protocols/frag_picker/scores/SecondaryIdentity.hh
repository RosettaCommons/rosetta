// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/SecondaryIdentity.hh
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_SecondaryIdentity_hh
#define INCLUDED_protocols_frag_picker_scores_SecondaryIdentity_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  SequenceIdentity score counts how many residues share the same secondary structure
/// @details Resulting score is the number of residues
/// on corresponding positions in a vall fragment and a query sequence
/// that are annotated with the same secondary structure
class SecondaryIdentity: public FragmentScoringMethod {
public:
	/// @brief will be added if the two residues are in the same secondary structure
	/// @details the REWARD value is 0 because we are minimizing the score and the perfect match should be 0
	static const int REWARD = 0;
	/// @brief will be added if the two residues are NOT in the same secondary structure
	static const int PENALTY = 1;

	SecondaryIdentity(Size priority, Real lowest_acceptable_value, bool use_lowest,
			std::string& query_secondary_structure,std::string prediction_name) :
		FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest,
				"SecondaryIdentity"), query_(query_secondary_structure),prediction_name_(prediction_name) {
	}

	/// @brief Computes the score
	virtual bool score(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief prints a detailed explanation how a fragment score has been computed
	/// @details besides extensive output, the method should return the same result as score()
	bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);

        inline std::string& get_prediction_name() { return prediction_name_; };
private:
	std::string& query_;
	std::string prediction_name_;
};

/// @brief  Maker class that produces a new SecondaryIdentity object
class MakeSecondaryIdentity: public MakeFragmentScoringMethod {
public:

	MakeSecondaryIdentity() :
		MakeFragmentScoringMethod("SecondaryIdentity") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value,
			bool use_lowest, FragmentPickerOP picker, std::string prediction_id) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new SecondaryIdentity(priority,
				lowest_acceptable_value, use_lowest, picker->get_query_ss_string(
						prediction_id),prediction_id) );
	}
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_SecondaryIdentity_HH */
