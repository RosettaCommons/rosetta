// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/Phi.hh
/// @brief  Object that scores a fragment by predicted Phi torsions
/// @author David E Kim

#ifndef INCLUDED_protocols_frag_picker_scores_Phi_hh
#define INCLUDED_protocols_frag_picker_scores_Phi_hh

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

using namespace core;

/// @brief  scores a fragment by its predicted phi similarity
class Phi: public CachingScoringMethod {
public:

	/// @brief  creates a predicted phi-based scoring function.
	Phi(Size priority, Real lowest_acceptable_value, bool use_lowest,
				std::string & fastaQuerySequence, utility::vector1<core::Real> & query_phi_prediction,
				utility::vector1<core::Real> & query_phi_prediction_conf) :
			CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
				"Phi"),  query_(fastaQuerySequence) {

				if (query_.length() != query_phi_prediction.size())
					utility_exit_with_message("Query length does not match predicted phi values");

				query_phi_prediction_ = query_phi_prediction;
				query_phi_prediction_conf_ = query_phi_prediction_conf;
	}

	~Phi() {};

	void do_caching(VallChunkOP);
	void clean_up() {};
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	std::string & query_;
  utility::vector1<Real> query_phi_prediction_;
  utility::vector1<Real> query_phi_prediction_conf_;
};

/// @brief  Maker class that produces a new Phi object
class MakePhi: public MakeFragmentScoringMethod {
public:

	MakePhi() :
		MakeFragmentScoringMethod("Phi") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value, bool use_lowest,
				FragmentPickerOP picker, std::string) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new Phi(priority,
        lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), picker->get_query_phi_prediction(),
				picker->get_query_phi_prediction_conf()) );
  }

};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_Phi_hh
