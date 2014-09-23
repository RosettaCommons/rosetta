// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite && is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/PValuedFragmentScoreManager.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// package headers
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
// AUTO-REMOVED #include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
// AUTO-REMOVED #include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/AdaptiveScoreHistogram.hh>
#include <protocols/frag_picker/scores/PValuedFragmentScoreManager.hh>

#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// AUTO-REMOVED #include <utility/io/izstream.hh>

// C++
#include <algorithm>
// AUTO-REMOVED #include <utility>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief print the scores and associated p-values
void PValuedFragmentScoreManager::describe_fragments(utility::vector1<std::pair<
		FragmentCandidateOP, scores::FragmentScoreMapOP> > const& pairs,
		std::ostream& out) {

	using namespace ObjexxFCL::format;

        out << "#" << RJ(10, "query_pos ");
        out << RJ(10, "vall_pos ");
	out << RJ(6, "pdbid");
        out << " c ss ";
	utility::vector1<Size> w(scores_.size()*2+2);
	for (Size i = 1; i <= scores_.size(); i++) {
		w[i*2 - 1] = width_.find(scores_[i])->second;
		out << " " << scores_[i]->get_score_name();
		if (scores_[i]->get_score_name().length() > w[i*2 - 1])
			w[i*2 - 1] = scores_[i]->get_score_name().length();
		w[i*2] = width_.find(scores_[i])->second;
		out << " p_" << scores_[i]->get_score_name();
		if (scores_[i]->get_score_name().length() + 2 > w[i*2])
			w[i*2] = scores_[i]->get_score_name().length() + 2;
	}
	w[scores_.size()*2+1] = 7; // 7 characters for the total score
	w[scores_.size()*2+2] = 7; // 7 characters for the total p-value
	out << "  TOTAL    p-TOTAL  FRAG_ID"<<std::endl;

	for (Size iF = 1; iF <= pairs.size(); ++iF) {

		FragmentCandidateOP fr = pairs[iF].first;
		FragmentScoreMapOP sc = pairs[iF].second;

		out << " " << I(10, fr->get_first_index_in_query());
		out << " " << I(10, fr->get_first_index_in_vall());
		out << " " << RJ(5, fr->get_pdb_id());
		out << " " << fr->get_chain_id();
		out << " " << fr->get_middle_ss();
		Real p_val_sum = 0.0;
		for (Size i = 1; i <= scores_.size(); i++) {
			Size p = precision_.find(scores_[i])->second;
			out << " " << F(w[i*2-1], p, sc->get_score_components()[i]);
			out << " " << F(w[i*2], 3, statistics_[fr->get_first_index_in_query()][i]->p_value(sc->get_score_components()[i]));
			p_val_sum += statistics_[fr->get_first_index_in_query()][i]->p_value(sc->get_score_components()[i]);
		}
		out << F(w[scores_.size()*2+1],TOTAL_PRECISION,total_score(sc));
		out << F(w[scores_.size()*2+2],TOTAL_PRECISION,p_val_sum);

    		out << I(10, fr->key() ) << std::endl;
	}
}


bool PValuedFragmentScoreManager::score_fragment(FragmentCandidateOP candidate,
		FragmentScoreMapOP empty_map) {

	FragmentScoreManager::score_fragment(candidate,empty_map);

	Size pos = candidate->get_first_index_in_query();
	if(statistics_.size() < pos) {
	    Size size_is = statistics_.size();
	    statistics_.resize( pos );
	    for(Size k=size_is+1;k<=pos;k++) {
		for(Size l=1;l<=scores_.size();l++)
		    statistics_[k].push_back(utility::pointer::shared_ptr<class protocols::frag_picker::scores::AdaptiveScoreHistogram>( new AdaptiveScoreHistogram(0.01,1) ));
	    }
	}
	for (Size iScore = 1; iScore <= empty_map->size(); iScore++) {
	    statistics_[pos][iScore]->insert( empty_map->at(iScore) );
	}

	return true;
}

bool PValuedFragmentScoreManager::score_fragment_from_cache(FragmentCandidateOP candidate,
		FragmentScoreMapOP empty_map) {

	FragmentScoreManager::score_fragment_from_cache(candidate,empty_map);

	Size pos = candidate->get_first_index_in_query();
	if(statistics_.size() < pos) {
	    Size size_is = statistics_.size();
	    statistics_.resize( pos );
	    for(Size k=size_is+1;k<=pos;k++) {
		for(Size l=1;l<=scores_.size();l++)
		    statistics_[k].push_back(utility::pointer::shared_ptr<class protocols::frag_picker::scores::AdaptiveScoreHistogram>( new AdaptiveScoreHistogram(0.01,1) ));
	    }
	}

	for (Size iScore = 1; iScore <= empty_map->size(); iScore++) {
	    statistics_[pos][iScore]->insert( empty_map->at(iScore) );
	}

	return true;
}


} // scores
} // frag_picker
} // protocols


