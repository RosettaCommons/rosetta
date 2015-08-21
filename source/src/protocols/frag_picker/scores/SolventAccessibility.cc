// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/SolventAccessibility.cc
/// @brief  a base class for fragment scoring
/// @author David E Kim

#include <protocols/frag_picker/scores/SolventAccessibility.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

void SolventAccessibility::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool SolventAccessibility::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	return score( fragment, scores );

}

bool SolventAccessibility::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {
	Real totalScore = 0;
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);
		// just the difference
		totalScore += fabs( predicted_sa_norm_[i + f->get_first_index_in_query() - 1] - r->sa_norm() );
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}


} // scores
} // frag_picker
} // protocols


