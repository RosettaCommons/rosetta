// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/HydrophobicityProfileSimilarity.hh
/// @brief  Object that scores a fragment by hydrophobicity
/// @author David E Kim

#ifndef INCLUDED_protocols_frag_picker_scores_HydrophobicityProfileSimilarity_hh
#define INCLUDED_protocols_frag_picker_scores_HydrophobicityProfileSimilarity_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// type headers
#include <core/types.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>

#include <map>

namespace protocols {
namespace frag_picker {
namespace scores {


/// @brief  scores a fragment by its hydrophobicity similarity
class HydrophobicityProfileSimilarity: public CachingScoringMethod {
public:

	/// @brief  creates a hydrophobicity-based scoring function.
	HydrophobicityProfileSimilarity(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		std::string & fastaQuerySequence, core::sequence::SequenceProfileOP query_profile) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
		"HydrophobicityProfileSimilarity"),  query_(fastaQuerySequence) {
		if ( query_profile->size() != query_profile->length() ) {
			utility_exit_with_message("HydrophobicityProfileSimilarity needs a valid sequence profile.");
		}

		query_profile_ = query_profile;

		// initialize hydrophobic aa map
		// V, I, L, F, Y, W, M
		if ( is_hydrophobic_.empty() ) {
			is_hydrophobic_['A'] = false;
			is_hydrophobic_['C'] = false;
			is_hydrophobic_['D'] = false;
			is_hydrophobic_['E'] = false;
			is_hydrophobic_['F'] = true;
			is_hydrophobic_['G'] = false;
			is_hydrophobic_['H'] = false;
			is_hydrophobic_['I'] = true;
			is_hydrophobic_['K'] = false;
			is_hydrophobic_['L'] = true;
			is_hydrophobic_['M'] = true;
			is_hydrophobic_['N'] = false;
			is_hydrophobic_['P'] = false;
			is_hydrophobic_['Q'] = false;
			is_hydrophobic_['R'] = false;
			is_hydrophobic_['S'] = false;
			is_hydrophobic_['T'] = false;
			is_hydrophobic_['V'] = true;
			is_hydrophobic_['W'] = true;
			is_hydrophobic_['Y'] = true;
			is_hydrophobic_['a'] = false;
			is_hydrophobic_['c'] = false;
			is_hydrophobic_['d'] = false;
			is_hydrophobic_['e'] = false;
			is_hydrophobic_['f'] = true;
			is_hydrophobic_['g'] = false;
			is_hydrophobic_['h'] = false;
			is_hydrophobic_['i'] = true;
			is_hydrophobic_['k'] = false;
			is_hydrophobic_['l'] = true;
			is_hydrophobic_['m'] = true;
			is_hydrophobic_['n'] = false;
			is_hydrophobic_['p'] = false;
			is_hydrophobic_['q'] = false;
			is_hydrophobic_['r'] = false;
			is_hydrophobic_['s'] = false;
			is_hydrophobic_['t'] = false;
			is_hydrophobic_['v'] = true;
			is_hydrophobic_['w'] = true;
			is_hydrophobic_['y'] = true;
		}

	}

	~HydrophobicityProfileSimilarity() {};

	void do_caching(VallChunkOP);
	void clean_up() {};
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	std::string & query_;
	core::sequence::SequenceProfileOP query_profile_;
	static std::map<char,bool> is_hydrophobic_;
};

/// @brief  Maker class that produces a new HydrophobicityProfileSimilarity object
class MakeHydrophobicityProfileSimilarity: public MakeFragmentScoringMethod {
public:

	MakeHydrophobicityProfileSimilarity() :
		MakeFragmentScoringMethod("HydrophobicityProfileSimilarity") {
	}

	FragmentScoringMethodOP make(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
		FragmentPickerOP picker, std::string) {
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new HydrophobicityProfileSimilarity(priority,
			lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), picker->get_query_seq()) );
	}

};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_HydrophobicityProfileSimilarity_hh
