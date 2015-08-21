// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/HydrophobicitySimilarity.cc
/// @brief  a base class for fragment scoring
/// @author David E Kim

#include <protocols/frag_picker/scores/HydrophobicitySimilarity.hh>

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

void HydrophobicitySimilarity::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool HydrophobicitySimilarity::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	return score( fragment, scores );

}

bool HydrophobicitySimilarity::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {
	Real totalScore = 0;
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);
		// scoring based on SEGMER hydrophobicity energy term
		if ( is_hydrophobic_[r->aa()] && is_hydrophobic_[query_[i + f->get_first_index_in_query() - 2]] ) {
			continue; // no penalty if both query and vall residue are hydrophobic
		} else if ( r->aa() == query_[i + f->get_first_index_in_query() - 2] ) {
			if ( r->aa() == 'G' || r->aa() == 'g' || r->aa() == 'P' || r->aa() == 'p' ) {
				continue; // no penalty for identical G or P
			} else {
				totalScore += 0.3; // small penalty if residues are equal
			}
		} else {
			totalScore += 1.0; // penalty if different and not hydrophobic
		}
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

std::map<char,bool> HydrophobicitySimilarity::is_hydrophobic_;

} // scores
} // frag_picker
} // protocols


