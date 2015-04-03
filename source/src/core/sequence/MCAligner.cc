// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MCAligner.cc
/// @brief class definition for a class that aligns two Sequence objects using
/// a stochastic version of the Needleman-Wunsch algorithm.
/// @author James Thompson

#include <core/types.hh>

#include <core/sequence/Aligner.hh>
#include <core/sequence/MCAligner.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <numeric/random/random.hh>

#include <utility/vector1.hh>

#include <iostream>
#include <string>

namespace core {
namespace sequence {


	SequenceAlignment MCAligner::align(
		SequenceOP seq_y,
		SequenceOP seq_x,
		ScoringSchemeOP ss
	) {
		using std::exp;
		// replace these with clones, fix up match/mismatch calculation
		SequenceOP new_seq_x = seq_x->clone();
		SequenceOP new_seq_y = seq_y->clone();

		new_seq_x->insert_gap(0);
		new_seq_y->insert_gap(0);

		//std::cout << "originals are: " << (*seq_x) << std::endl << (*seq_y) << std::endl;
		//std::cout << "copies are: " << (*new_seq_x) << std::endl << (*new_seq_y) << std::endl;

		DP_Matrix scores( new_seq_x->length(), new_seq_y->length() );
		//scores.xlab( new_seq_x->sequence_vec() );
		//scores.ylab( new_seq_y->sequence_vec() );

		// std::cout << scores << std::endl;

		// matrix initialization
		NWAligner::init_matrix( new_seq_y->length(), new_seq_x->length(), ss, scores );
		// std::cout << seq_x << std::endl << seq_y << std::endl;
		// std::cout << scores << std::endl;

		CellOP best_cell( new Cell( 0.0 ) );
		for ( Size y = 2; y <= new_seq_y->length(); ++y ) {
			for ( Size x = 2; x <= new_seq_x->length(); ++x ) {

				// score this position
				Real l_gap_penalty( ss->gap_open() ), u_gap_penalty( ss->gap_open() );
				if ( scores(x,y-1)->came_from() == above ) u_gap_penalty = ss->gap_extend();
				if ( scores(x-1,y)->came_from() == left  ) l_gap_penalty = ss->gap_extend();

				Real single_mm    = ss->score( new_seq_x, new_seq_y, x, y );

				Real u_gap = scores( x, y-1   )->score() + u_gap_penalty;
				Real l_gap = scores( x-1,   y )->score() + l_gap_penalty;
				Real mm    = scores( x-1, y-1 )->score() + ss->score( new_seq_x, new_seq_y, x, y );

				// convert u_gap, l_gap and mm to probs using boltzmann averaging (based on kT)
				Real u_gap_prob( exp( -1 * u_gap_penalty / kT() ) );
				Real l_gap_prob( exp( -1 * l_gap_penalty / kT() ) );
				Real mm_prob   ( exp( -1 * single_mm     / kT() ) );

				Real const partition( u_gap_prob + l_gap_prob + mm_prob );
				Real const rand_prob( numeric::random::rg().uniform() );
				//Real const rand_prob( numeric::random::uniform() );
				u_gap_prob /= partition;
				l_gap_prob /= partition;
				mm_prob    /= partition;

				// think of rand_prob as a dart that hits one of three places:
				// (0,u_gap_prob)          => align with u_gap
				// (u_gap_prob,l_gap_prob) => align with l_gap
				// (l_gap_prob,1)          => align with mm

				CellOP current_cell = scores(x,y);
				if ( rand_prob < u_gap_prob ) {
					// std::cout << "came from diagonal with a score of " << mm << std::endl;
					current_cell->score( mm );
					current_cell->next( scores(x-1,y-1) );
					current_cell->came_from( diagonal );
				} else if ( rand_prob < l_gap_prob ) { // && rand_prob > 0, handled earlier
					// std::cout << "came from left with a score of " << l_gap << std::endl;
					current_cell->score( l_gap );
					current_cell->next( scores(x-1,y) );
					current_cell->came_from( left );
				} else { // if rand_prob > l_gap_prob && rand_prob < 1
					// std::cout << "came from above with a score of " << u_gap << std::endl;
					current_cell->score( u_gap );
					current_cell->next( scores(x,y-1) );
					current_cell->came_from( above );
				}
				// if ( current_cell->score() > best_cell->score() ) best_cell = current_cell;
			} // x
		} // y

		best_cell = scores( scores.rows(),scores.cols() ); // start in lower right-hand corner for NW
		SequenceAlignment test_alignment = traceback(
			new_seq_x,
			new_seq_y,
			scores,
			best_cell
		);

		scores.clear();

		return test_alignment;
	}  // align

} // sequence
} // core
