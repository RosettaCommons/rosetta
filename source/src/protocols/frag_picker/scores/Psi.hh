// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/Psi.hh
/// @brief  Object that scores a fragment by predicted Psi torsions
/// @author David E Kim

#ifndef INCLUDED_protocols_frag_picker_scores_Psi_hh
#define INCLUDED_protocols_frag_picker_scores_Psi_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// utility headers
#include <utility/exit.hh>

// type headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

/// @brief  scores a fragment by its predicted psi similarity
class Psi: public CachingScoringMethod {
public:

	/// @brief  creates a predicted psi-based scoring function.
	Psi(Size priority, Real lowest_acceptable_value, bool use_lowest,
		std::string & fastaQuerySequence, utility::vector1<core::Real> & query_psi_prediction,
		utility::vector1<core::Real> & query_psi_prediction_conf) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
		"Psi"),  query_(fastaQuerySequence) {

		if ( query_.length() != query_psi_prediction.size() ) {
			utility_exit_with_message("Query length does not match predicted psi values");
		}

		query_len_ = query_.length();
		query_psi_prediction_ = query_psi_prediction;
		query_psi_prediction_conf_ = query_psi_prediction_conf;
	}

	~Psi() {};

	void do_caching(VallChunkOP);
	void clean_up() {};
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	std::string & query_;
	Size query_len_;
	utility::vector1<Real> query_psi_prediction_;
	utility::vector1<Real> query_psi_prediction_conf_;
};

/// @brief  Maker class that produces a new Psi object
class MakePsi: public MakeFragmentScoringMethod {
public:

	MakePsi() :
		MakeFragmentScoringMethod("Psi") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value, bool use_lowest,
		FragmentPickerOP picker, std::string) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new Psi(priority,
			lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), picker->get_query_psi_prediction(),
			picker->get_query_psi_prediction_conf()) );
	}

};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_Psi_hh
