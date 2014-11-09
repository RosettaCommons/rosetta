// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite && is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.cc
/// @brief  scores a fragment by an amino acid sequence identity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScore.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>
#include <basic/Tracer.hh>

// option key includes
// AUTO-REMOVED #include <core/init/init.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// utils
// AUTO-REMOVED #include <basic/prof.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trProfScore(
		"protocols.frag_picker.scores.ProfileScore");

ProfileScore::~ProfileScore() {}

void ProfileScore::do_caching(VallChunkOP chunk) {

	std::string tmp = chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;
	Size size_q = query_profile_->length();
	core::sequence::SequenceProfileOP chunk_profile = chunk->get_profile();
	assert( chunk_profile->length() !=0 );
	trProfScore.Debug << "caching profile score for " << chunk->get_pdb_id()
			<< " of size " << chunk->size() << std::endl;
	assert( chunk->size() == chunk_profile->length() );
	for (Size i = 1; i <= size_q; ++i) {
		for (Size j = 1; j <= chunk->size(); ++j) {
			scores_[i][j] = profile_scoring_->score(query_profile_,
					chunk_profile, i, j);
		}
	}

	trProfScore.Debug << "precomputed matrix of scores " << scores_.size()
			<< "x" << chunk->size() << std::endl;
}

bool ProfileScore::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

/*	std::string tmp = f->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(f->get_chunk());
*/

	Real totalScore = 0;
	for (Size i = 1; i <= f->get_length(); i++) {
		assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		assert(f->get_first_index_in_vall()
				+ i - 1<= scores_[1].size());
		totalScore
				+= scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall()
						+ i - 1];
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

bool ProfileScore::describe_score(FragmentCandidateOP f,
		FragmentScoreMapOP empty_map, std::ostream& out) {

	Real totalScore = 0;
	core::sequence::SequenceProfileOP chunk_profile =
			f->get_chunk()->get_profile();
	Size firstQ = f->get_first_index_in_query();
	Size firstV = f->get_first_index_in_vall();

/*
#ifndef WIN32
	out << "\nvall: " << f->get_chunk()->get_pdb_id() << " " << I(5,
			f->get_first_index_in_vall()) << " " << out << "\nquery " << I(5,
			f->get_first_index_in_query()) << "\n";
#else */
	// In VS2013, can't do: out << out
	// libc++ with c++11 can't do this either. I have commented the non-working version of this and contacted Dominik about this.
	out << "\nvall: " << f->get_chunk()->get_pdb_id() << " " << I(5,
			f->get_first_index_in_vall()) << "\nquery " << I(5,
			f->get_first_index_in_query()) << "\n";
// #endif

	for (Size i = 1; i <= f->get_length(); ++i) {
		out << "V row: " << F(5, 3, chunk_profile->prof_row(firstV + i - 1)[1]);
		for (Size j = 2; j <= 20; ++j) {
			out << " " << F(5, 3, chunk_profile->prof_row(firstV + i - 1)[j]);
		}
		out << "\n";
		out << "Q row: "
				<< F(5, 3, query_profile_->prof_row(firstQ + i - 1)[1]);
		for (Size j = 2; j <= 20; ++j) {
			out << " " << F(5, 3, query_profile_->prof_row(firstQ + i - 1)[j]);
		}
		out << "\n";
		out << "Position score: " << scores_[firstQ + i - 1][firstV + i - 1]
				<< "\n";
		totalScore
				+= scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall()
						+ i - 1];
	}
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		out << "Total score " << F(5, 3, totalScore) << " ACCEPTED"
				<< std::endl;
	else
		out << "Total score " << F(5, 3, totalScore) << " REJECTED"
				<< std::endl;
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

FragmentScoringMethodOP MakeProfileScore::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	if (option[frags::scoring::profile_score].user()) {
		Size len = picker->get_vall()->get_largest_chunk_size();
		core::sequence::ScoringSchemeFactory ssf;
		core::sequence::ScoringSchemeOP ss(ssf.get_scoring_scheme(
				option[frags::scoring::profile_score]()));
		trProfScore << "Profile scoring method is: "
				<< option[frags::scoring::profile_score]() << std::endl;
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new ProfileScore(priority,
				lowest_acceptable_value, use_lowest, picker->get_query_seq(), ss, len) );
	}
	utility_exit_with_message(
			"[ERROR] Undefined profile scoring method. Provide it with frags::scoring::profile_score flag");

	return NULL;
}

} //scores
} // frag_picker
} // protocols
