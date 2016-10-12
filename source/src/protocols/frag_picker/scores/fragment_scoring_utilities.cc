// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/fragment_scoring_utilities.hh
/// @brief  functions and data types common for fragment scoring
/// @author Dominik Gront

#include <protocols/frag_picker/scores/fragment_scoring_utilities.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

void create_cache(utility::vector1<core::Size> & frag_sizes,core::Size query_len,core::Size longest_vall_chunk,utility::vector1<Matrix> & cache) {

	std::sort(frag_sizes.begin(),frag_sizes.end());
	cache.resize(frag_sizes[frag_sizes.size()]);
	for ( core::Size i=1; i<=frag_sizes.size(); i++ ) {
		Matrix m;
		allocate_matrix(query_len,longest_vall_chunk,m);
		cache[frag_sizes[i]] = m;
	}
}

void allocate_matrix(core::Size i_size,core::Size j_size,Matrix & dst) {

	dst.clear();
	for ( core::Size i = 1; i <= i_size; i++ ) {
		utility::vector1<core::Real> row(j_size);
		dst.push_back( row );
	}
}

void do_one_line(core::Size start_i,core::Size start_j,Matrix & small_scores,core::Size frag_len,Matrix & frag_scores) {

	core::Size stop_i = start_i + frag_len - 1;
	core::Size stop_j = start_j + frag_len - 1;
	core::Real last_score = small_scores[start_i][start_j];

	for ( core::Size i=1; i<frag_len; i++ ) {
		last_score += small_scores[start_i+i][start_j+i];
	}
	frag_scores[start_i][start_j] = last_score;
	int max_steps = std::min((int) small_scores.size() - (int) start_i - (int) frag_len,
		(int) small_scores[1].size() - (int) start_j - (int) frag_len) + 1;
	int cnt = 1;
	while ( cnt<=max_steps ) {
		stop_i++;
		stop_j++;
		last_score += small_scores[stop_i][stop_j] - small_scores[start_i][start_j];
		start_i++;
		start_j++;
		frag_scores[start_i][start_j] = last_score;
		cnt++;
	}
}

void rolling_score(Matrix & small_scores,core::Size frag_len,Matrix & frag_scores) {

	do_one_line(1,1,small_scores,frag_len,frag_scores);
	for ( core::Size i=2; i<=small_scores.size()-frag_len+1; i++ ) {
		do_one_line(i,1,small_scores,frag_len,frag_scores);
	}
	for ( core::Size i=2; i<=small_scores[1].size()-frag_len+1; i++ ) {
		do_one_line(1,i,small_scores,frag_len,frag_scores);
	}
}


} // scores
} // frag_picker
} // protocols

