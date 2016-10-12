// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/RDCScore.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Ray Wang ( wangy@u.washington.edu )

#ifndef INCLUDED_protocols_frag_picker_scores_RDCScore_hh
#define INCLUDED_protocols_frag_picker_scores_RDCScore_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/scoring/ResidualDipolarCoupling.hh>

// type headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief  scores a fragment by its crmsd to the given reference structure
class RDCScore: public CachingScoringMethod {
public:
	/// @brief  creates a RDC-based scoring function.
	RDCScore( core::Size, core::Real, bool);

	void do_caching( VallChunkOP );
	void clean_up();
	bool cached_score( FragmentCandidateOP, FragmentScoreMapOP );
	bool score(FragmentCandidateOP, FragmentScoreMapOP);


private:
	core::scoring::ResidualDipolarCouplingOP rdc_file_;
	core::scoring::ResidualDipolarCoupling::RDC_lines rdc_raw_data_;
};

/// @brief  Maker class that produces a new RDCScore object
class MakeRDCScore: public MakeFragmentScoringMethod {
public:

	MakeRDCScore() :
		MakeFragmentScoringMethod("RDCScore") {
	}

	FragmentScoringMethodOP make( core::Size, core::Real, bool, FragmentPickerOP, std::string );
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_RDCScore_HH
