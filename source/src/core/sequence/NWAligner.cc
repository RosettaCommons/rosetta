// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NWAligner.cc
/// @brief class definition for a class that aligns two Sequence objects using
/// the Needleman-Wunsch alignment algorithm.
/// @author James Thompson

#include <core/types.hh>

#include <core/sequence/Aligner.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

void NWAligner::init_matrix(
	Size y_len,
	Size x_len,
	ScoringSchemeOP ss,
	DP_Matrix & scores
) {
	// matrix initialization
	for ( Size y = 1; y <= y_len; ++y ) {
		if ( y == 1 ) {
			scores(1,y)->score( 0.0 );
			scores(1,y)->came_from( end );
		} else {
			scores(1,y)->score( ss->gap_open() + (y-2) * ss->gap_extend() );
			scores(1,y)->came_from( end );
			scores(1,y)->next( scores(1,y-1) );
		}
	} // for y

	for ( Size x = 1; x <= x_len; ++x ) {
		if ( x == 1 ) {
			scores(x,1)->score( 0.0 );
			scores(x,1)->came_from( end );
		} else {
			scores(x,1)->score( ss->gap_open() + (x-2) * ss->gap_extend() );
			scores(x,1)->came_from( end );
			scores(x,1)->next( scores(x-1,1) );
		}
	} // for x
} // init_matrix

SequenceAlignment NWAligner::align(
	SequenceOP seq_y,
	SequenceOP seq_x,
	ScoringSchemeOP ss
) {
	validate_input( seq_y, seq_x, ss );

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
	init_matrix( new_seq_y->length(), new_seq_x->length(), ss, scores );

	// std::cout << seq_x << std::endl << seq_y << std::endl;
	// std::cout << scores << std::endl;

	CellOP best_cell( new Cell( 0.0 ) );
	for ( Size y = 2; y <= new_seq_y->length(); ++y ) {
		for ( Size x = 2; x <= new_seq_x->length(); ++x ) {

			// score this position
			Real l_gap_penalty( ss->gap_open() ), u_gap_penalty( ss->gap_open() );
			if ( scores(x,y-1)->came_from() == above ) u_gap_penalty = ss->gap_extend();
			if ( scores(x-1,y)->came_from() == left  ) l_gap_penalty = ss->gap_extend();

			Real u_gap = scores( x, y-1   )->score() + u_gap_penalty;
			Real l_gap = scores( x-1,   y )->score() + l_gap_penalty;
			Real mm    = scores( x-1, y-1 )->score() + ss->score( new_seq_x, new_seq_y, x, y );

			// set up pointers properly
			CellOP current_cell = scores(x,y);
			if ( mm >= l_gap && mm >= u_gap ) {
				// std::cout << "came from diagonal with a score of " << mm << std::endl;
				current_cell->score( mm );
				current_cell->next( scores(x-1,y-1) );
				current_cell->came_from( diagonal );
			} else if ( l_gap >= mm && l_gap >= u_gap ) {
				// std::cout << "came from left with a score of " << l_gap << std::endl;
				current_cell->score( l_gap );
				current_cell->next( scores(x-1,y) );
				current_cell->came_from( left );
			} else { // if ( u_gap >= mm && u_gap >= l_gap ) {
				// std::cout << "came from above with a score of " << u_gap << std::endl;
				current_cell->score( u_gap );
				current_cell->next( scores(x,y-1) );
				current_cell->came_from( above );
			}
			// if ( current_cell->score() > best_cell->score() ) best_cell = current_cell;
		} // x
	} // y

	best_cell = scores( scores.rows(),scores.cols() ); // start in lower right-hand corner for NW

	SequenceAlignment aln = traceback(
		new_seq_x,
		new_seq_y,
		scores,
		best_cell
	);

	scores.clear();

	return aln;
}  // align

} // sequence
} // core
