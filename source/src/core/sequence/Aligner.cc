// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Aligner.cc
/// @brief class definition for a class that aligns two Sequence objects using
/// dynamic programming algorithms.
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Aligner.hh>

#include <core/sequence/DP_Matrix.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>


#include <string>

#include <utility/exit.hh>

#include <utility/vector1.hh>

#include <core/sequence/ScoringScheme.fwd.hh>

#ifdef WIN32
#include <core/sequence/ScoringScheme.hh>
#endif


namespace core {
namespace sequence {

void Aligner::validate_input(
	SequenceOP seq_y,
	SequenceOP seq_x,
	ScoringSchemeOP ss
) {
	runtime_assert( ss    != nullptr );
	runtime_assert( seq_x != nullptr );
	runtime_assert( seq_y != nullptr );

	// must have ungapped sequences!
	runtime_assert( seq_y->ungapped_length() == seq_y->length() );
	runtime_assert( seq_x->ungapped_length() == seq_x->length() );
}

SequenceAlignment Aligner::traceback(
	SequenceOP seq_x,
	SequenceOP seq_y,
	DP_Matrix /*matrix*/,
	CellOP start
) {
	// traceback
	CellOP current_cell = start;
	std::string aligned_seq_x(""), aligned_seq_y("");

	while ( true ) {
		Size current_x = current_cell->x();
		Size current_y = current_cell->y();

		//std::cout << "at " << current_x << "," << current_y << std::endl;

		if ( current_cell->came_from() == diagonal ) {
			//std::cout << " came from diagonal from score of " << current_cell->next()->score() << std::endl;
			aligned_seq_x = (*seq_x)[ current_x ] + aligned_seq_x;
			aligned_seq_y = (*seq_y)[ current_y ] + aligned_seq_y;
		} else if ( current_cell->came_from() == left ) {
			//std::cout << " came from left from score of " << current_cell->next()->score() << std::endl;
			aligned_seq_x = (*seq_x)[current_x] + aligned_seq_x;
			aligned_seq_y = '-' + aligned_seq_y;
			seq_y->insert_gap( current_y + 1 );
			// std::cout << "seq_x[" << current_x << "] = " << seq_x[current_x] << std::endl;
		} else if ( current_cell->came_from() == above ) {
			//std::cout << " came from above from score of " << current_cell->next()->score() << std::endl;
			aligned_seq_x = '-' + aligned_seq_x;
			aligned_seq_y = (*seq_y)[current_y] + aligned_seq_y;
			seq_x->insert_gap( current_x + 1 );
		} else {
			//std::string const msg(
			// "Unhandled case in traceback, not pointing to anything (" +
			// "\n"
			// //string_of(current_x) + "," +
			// //string_of(current_y) + ")!\n"
			//);
			//utility_exit_with_message( msg );
			utility_exit_with_message( "Error in traceback: pointer doesn't go anywhere!\n" );
		}

		if ( current_cell->next()->came_from() == end ) {
			break;
		}

		current_cell = current_cell->next();
		//std::cout << aligned_seq_x << std::endl << aligned_seq_y << std::endl << std::endl;
	} // while ( current_cell->next() != 0 )

	//std::cout << std::endl << (*seq_x) << std::endl << (*seq_y) << std::endl;
	//std::cout << matrix << std::endl;

	SequenceAlignment alignment;
	// set starting point for both sequences. Don't forget to add whatever offset existed
	// in Sequence.start() coming into this function. Also, subtract one for the extra gap
	// inserted at the beginning of new_seq_x and new_seq_y.
	SequenceOP seq_x_clone = seq_x->clone();
	SequenceOP seq_y_clone = seq_y->clone();

	seq_x_clone->start( current_cell->x() - 2 + seq_x->start() );
	seq_y_clone->start( current_cell->y() - 2 + seq_y->start() );

	seq_x_clone->sequence( aligned_seq_x );
	seq_y_clone->sequence( aligned_seq_y );

	alignment.add_sequence( seq_y_clone );
	alignment.add_sequence( seq_x_clone );

	alignment.remove_gapped_positions();
	alignment.score( start->score() );

	return alignment;
} // traceback

} // sequence
} // core
