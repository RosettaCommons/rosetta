// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DPScoringScheme.cc
/// @brief method implementations for DPScoringScheme class.
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/DPScoringScheme.hh>

#include <utility/exit.hh>
#include <string>


#include <utility/vector1.hh>


namespace core {
namespace sequence {

Real DPScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	Size pos1,
	Size pos2
) {
	SequenceProfileOP prof1 = SequenceProfileOP(
		utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq1 )
	);
	SequenceProfileOP prof2 = SequenceProfileOP(
		utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq2 )
	);

	runtime_assert( pos1 <= prof1->length() );
	runtime_assert( pos2 <= prof2->length() );
	runtime_assert( prof1->prof_row(pos1).size() == prof2->prof_row(pos2).size() );

	// compare the two profiles using a dot-product between the two profile
	// vectors.
	Size n_aa( prof1->prof_row(pos1).size() );
	Real score( 0.0 );
	for ( Size i = 1; i <= n_aa; ++i ) {
		Real const & p1_num( prof1->prof_row(pos1)[i] );
		Real const & p2_num( prof2->prof_row(pos2)[i] );
		if ( is_good( p1_num ) && is_good( p2_num ) && i == 2 ) {
		//if ( is_good( p1_num ) && is_good( p2_num ) ) {
			//for ( Size jj = 1; jj < i; ++jj ) std::cout << " ";
			//std::cout << "comparing " << p1_num << " and " << p2_num << std::endl;
			//score += std::min( std::abs( p1_num * p2_num ), 1.0 );
			score += p1_num * p2_num;
		}
	}

	//score = score - 6;
	//std::cout << "score(" << pos1 << "," << pos2 << ") = " << score << std::endl << std::endl;

	return score;
} // score

} // sequence
} // core
