// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/SolventAccessibility.hh
/// @brief  Object that scores a fragment by predicted solvent accessibility
/// @author David E Kim

#ifndef INCLUDED_protocols_frag_picker_scores_SolventAccessibility_hh
#define INCLUDED_protocols_frag_picker_scores_SolventAccessibility_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/Faraggi_SA.hh>

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

/// @brief  scores a fragment by its predicted solvent accessibility
class SolventAccessibility: public CachingScoringMethod {
public:

	/// @brief  creates a predicted solvent accessibility based scoring function.
	SolventAccessibility(Size priority, Real lowest_acceptable_value, bool use_lowest,
		std::string & fastaQuerySequence, utility::vector1<core::Real> & predicted_sa) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
		"SolventAccessibility"),  query_(fastaQuerySequence) {

		if ( query_.length() != predicted_sa.size() ) {
			utility_exit_with_message("Query length does not match predicted solvent accessiblity");
		}

		// get normalized ASA values
		//   just divide by the maximum values from Faraggi et al. Proteins 2008 (Table II)
		for ( Size i=1; i<=predicted_sa.size(); ++i ) {
			predicted_sa_norm_.push_back( predicted_sa[i]/protocols::frag_picker::sa_faraggi_max( query_[i-1] ) );
		}
	}

	~SolventAccessibility() {};

	void do_caching(VallChunkOP);
	void clean_up() {};
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	std::string & query_;
	utility::vector1<Real> predicted_sa_norm_;
};

/// @brief  Maker class that produces a new SolventAccessibility object
class MakeSolventAccessibility: public MakeFragmentScoringMethod {
public:

	MakeSolventAccessibility() :
		MakeFragmentScoringMethod("SolventAccessibility") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value, bool use_lowest,
		FragmentPickerOP picker, std::string) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new SolventAccessibility(priority,
			lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), picker->get_query_sa_prediction()) );
	}

};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_SolventAccessibility_hh
