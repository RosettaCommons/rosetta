// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SWAligner.hh
/// @brief class definition for a class that aligns two Sequence objects using
/// a ScoringScheme object and the Smith-Waterman alignment algorithm.
/// @author James Thompson

#include <core/types.hh>

#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

SequenceAlignment SWAligner::align(
	SequenceOP seq_y,
	SequenceOP seq_x,
	ScoringSchemeOP ss
) {
	validate_input( seq_y, seq_x, ss );

	SequenceOP new_seq_x = seq_x->clone();
	SequenceOP new_seq_y = seq_y->clone();

	new_seq_x->insert_gap(0);
	new_seq_y->insert_gap(0);

	DP_Matrix scores( new_seq_x->length(), new_seq_y->length() );
	//scores.xlab( new_seq_x->sequence_vec() );
	//scores.ylab( new_seq_y->sequence_vec() );

	// bail out of alignment if our score doesn't beat threshold
	Real const threshold( 0.0 );

	CellOP best_cell( new Cell( 0.0 ) );
	for ( Size y = 2; y <= new_seq_y->length(); ++y ) {
		for ( Size x = 2; x <= new_seq_x->length(); ++x ) {

			// score this position
			Real l_gap_penalty( ss->gap_open() ), u_gap_penalty( ss->gap_open() );
			if ( scores(x,y-1)->came_from() == above ) u_gap_penalty = ss->gap_extend();
			if ( scores(x-1,y)->came_from() == left  ) l_gap_penalty = ss->gap_extend();

			Real u_gap = scores( x,   y-1 )->score() + u_gap_penalty;
			Real l_gap = scores( x-1,   y )->score() + l_gap_penalty;
			Real mm    = scores( x-1, y-1 )->score() + ss->score( new_seq_x, new_seq_y, x, y );
			//Real mm    = scores( x-1, y-1 )->score() + ss->score( new_seq_x, new_seq_y, x, y ) - 0.45; // hacky

			//std::cout << mm << " = " << scores( x-1, y-1 )->score() << " + " << ss->score( new_seq_x, new_seq_y, x, y ) << std::endl;

			//std::cout << scores << std::endl;
			//std::cout << "(" << x << "," << y << ")"
			//     << " can choose from " << mm << "," << l_gap << "," << u_gap
			//     << std::endl;
			//std::cout << "mm = " << mm << "("
			// << scores( x-1, y-1 )->score() << "+"
			// << ss->score( new_seq_x, new_seq_y, x, y ) << ")"
			// << std::endl;
			CellOP current_cell = scores(x,y);
			if ( mm > l_gap && mm > u_gap && mm >= threshold ) {
				//if ( mm >= l_gap && mm >= u_gap && mm >= threshold ) {
				//std::cout << "came from diagonal with a score of " << mm << std::endl;
				current_cell->score( mm );
				current_cell->next( scores(x-1,y-1) );
				current_cell->came_from( diagonal );
			} else if ( l_gap >= mm && l_gap >= u_gap && l_gap >= threshold ) {
				//std::cout << "came from left with a score of " << l_gap << std::endl;
				current_cell->score( l_gap );
				current_cell->next( scores(x-1,y) );
				current_cell->came_from( left );
			} else if ( u_gap >= mm && u_gap >= l_gap && u_gap >= threshold ) {
				//std::cout << "came from above with a score of " << u_gap << std::endl;
				current_cell->score( u_gap );
				current_cell->next( scores(x,y-1) );
				current_cell->came_from( above );
			} else {
				current_cell->score( threshold );
				current_cell->came_from( end );
			}

			if ( current_cell->score() > best_cell->score() ) {
				best_cell = current_cell;
			}
		} // x
	} // y
	//std::cout << scores << std::endl;

	// traceback
	CellOP current_cell = best_cell;
	//std::cout << "best_cell = " << best_cell->x() << "," << best_cell->y() << std::endl;
	SequenceAlignment alignment = traceback(
		new_seq_x,
		new_seq_y,
		scores,
		best_cell
	);

	scores.clear();
	return alignment;
}  // align

} // sequence
} // core
