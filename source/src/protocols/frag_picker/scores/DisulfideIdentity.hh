// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/DisulfideIdentity.hh
/// @brief  a base class for fragment scoring
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_DisulfideIdentity_hh
#define INCLUDED_protocols_frag_picker_scores_DisulfideIdentity_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  DisulfideIdentity score counts identical residues
/// @details Resulting score is the number of identical residues
/// on corresponding positions in a vall fragment and a query sequence
class DisulfideIdentity: public FragmentScoringMethod {
public:

	DisulfideIdentity(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		std::string & fastaQueryDisulfide) :
		FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest,
		"DisulfideIdentity"), query_(fastaQueryDisulfide) {
	}

	bool score(FragmentCandidateOP f, FragmentScoreMapOP empty_map);

	/// @brief prints a detailed explanation how a fragment score has been computed
	/// @details besides extensive output, the method should return the same result as score()
	bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map,
		std::ostream& out);

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real minScoreAllowed_;
	std::string& query_;
};

/// @brief  Maker class that produces a new DisulfideIdentity object
class MakeDisulfideIdentity: public MakeFragmentScoringMethod {
public:

	MakeDisulfideIdentity() :
		MakeFragmentScoringMethod("DisulfideIdentity") {
	}

	FragmentScoringMethodOP make(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		FragmentPickerOP picker, std::string) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new DisulfideIdentity(priority,
			lowest_acceptable_value, use_lowest, picker->get_query_seq_string()) );
	}
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_DisulfideIdentity_HH */
