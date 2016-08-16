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

#ifndef INCLUDED_protocols_frag_picker_scores_fragment_scoring_utilities_hh
#define INCLUDED_protocols_frag_picker_scores_fragment_scoring_utilities_hh

// type headers
#include <core/types.hh>
#include <utility/vector1.hh>

#include <algorithm>

namespace protocols {
namespace frag_picker {
namespace scores {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

typedef utility::vector1< utility::vector1 < Real > > Matrix;


void do_one_line(Size start_i,Size start_j,Matrix & small_scores,Size frag_len,Matrix & frag_scores);
void rolling_score(Matrix & small_scores,Size frag_len,Matrix & frag_scores);
void create_cache(utility::vector1<Size> & frag_sizes,Size query_len,Size longest_vall_chunk,utility::vector1<Matrix> & cache);
void allocate_matrix(Size i_size,Size j_size,Matrix & dst);


} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_fragment_scoring_utilities_HH */
