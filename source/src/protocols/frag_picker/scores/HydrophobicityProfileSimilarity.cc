// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/HydrophobicityProfileSimilarity.cc
/// @brief  a base class for fragment scoring
/// @author David E Kim

#include <protocols/frag_picker/scores/HydrophobicityProfileSimilarity.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// mini headers
#include <core/chemical/AA.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

void HydrophobicityProfileSimilarity::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool HydrophobicityProfileSimilarity::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	return score( fragment, scores );

}

bool HydrophobicityProfileSimilarity::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	static const core::Real MIN_HYDROPHOBIC_PROBABILITY(0.75);

	core::Real totalScore = 0;
	VallChunkOP chunk = f->get_chunk();
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		core::Size qindex = i + f->get_first_index_in_query() - 1;
		if ( is_hydrophobic_[query_[qindex - 1]] ) {
			utility::vector1<core::Real> query_prof_row = query_profile_->prof_row(qindex);
			core::Real query_hydrophobic_sum = query_prof_row[core::chemical::aa_from_oneletter_code( 'F' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'I' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'L' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'M' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'V' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'W' )] +
				query_prof_row[core::chemical::aa_from_oneletter_code( 'Y' )];
			if ( query_hydrophobic_sum > MIN_HYDROPHOBIC_PROBABILITY ) {
				core::Size vallindex = f->get_first_index_in_vall() + i - 1;
				utility::vector1<core::Real> tmplt_prof_row = chunk->at(vallindex)->profile();
				core::Real tmplt_residue_hydrophobic_sum = tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'F' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'I' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'L' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'M' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'V' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'W' )] +
					tmplt_prof_row[core::chemical::aa_from_oneletter_code( 'Y' )];
				utility::vector1<core::Real> tmplt_prof_struct_row = chunk->at(vallindex)->profile_struct();
				tmplt_residue_hydrophobic_sum += tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'F' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'I' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'L' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'M' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'V' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'W' )] +
					tmplt_prof_struct_row[core::chemical::aa_from_oneletter_code( 'Y' )];
				tmplt_residue_hydrophobic_sum /= (core::Real) 2.0;
				if ( tmplt_residue_hydrophobic_sum < MIN_HYDROPHOBIC_PROBABILITY - 0.10 ) { // .1 buffer
					totalScore += 0.5; // penalty if avg sum of hydrophobic residue profile values from vall sequence and structure profiles less than query avg sum
				}
			}
		}
	}
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

std::map<char,bool> HydrophobicityProfileSimilarity::is_hydrophobic_;

} // scores
} // frag_picker
} // protocols


