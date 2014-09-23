// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ChemicalShiftScoringScheme.cc
/// @brief method implementations for ChemicalShiftScoringScheme class.
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ChemicalShiftSequence.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ChemicalShiftScoringScheme.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <iostream>
#include <string>

#include <complex>

#include <utility/vector1.hh>



namespace core {
namespace sequence {

Real ChemicalShiftScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	core::Size pos1,
	core::Size pos2
) {

	using core::Size;
	using core::Real;

	ChemicalShiftSequenceOP prof1 = ChemicalShiftSequenceOP(
		utility::pointer::dynamic_pointer_cast< core::sequence::ChemicalShiftSequence > ( seq1 )
	);
	ChemicalShiftSequenceOP prof2 = ChemicalShiftSequenceOP(
		utility::pointer::dynamic_pointer_cast< core::sequence::ChemicalShiftSequence > ( seq2 )
	);

	runtime_assert( pos1 <= prof1->length() );
	runtime_assert( pos2 <= prof2->length() );
	runtime_assert( prof1->prof_row(pos1).size() == prof2->prof_row(pos2).size() );

	//CA.q.dat.dev 0.828795509906
	//CB.q.dat.dev 0.969003290496
	//C.q.dat.dev  0.948093982461
	//HA.q.dat.dev 0.231641960268
	//HN.q.dat.dev 0.725943902221
	//N.q.dat.dev  2.78164308475

	Real const a(2), b(4); // constants for sigmoid smoothing trick

	Size n_aa( prof1->prof_row(pos1).size() );
	Size count( 0 );
	Real score( 0.0 );
	for ( Size i = 1; i <= n_aa; ++i ) {
		Real const & q_shift( prof1->prof_row(pos1)[i] );
		Real const & v_shift( prof2->prof_row(pos2)[i] );
		Real const & v_sigma( prof2->sigma(pos2,i) );

		if ( is_good( q_shift ) && is_good( v_shift ) ) {
			Real sig_diff(std::abs((q_shift - v_shift) / v_sigma ));
			Real sigmoid_diff( 1 / ( 1 + exp((a*sig_diff)-b) ) );

			//for ( Size jj = 1; jj < i; ++jj ) std::cout << " ";
			//std::cout << "comparing " << q_shift << " and " << v_shift <<
			//	" (sig_diff = " << sig_diff << ")" << std::endl <<
			//	" (v_sigma = " << v_sigma << ")" << std::endl;
			//for ( Size jj = 1; jj < i; ++jj ) std::cout << " ";
			//std::cout << "sigmoid_diff = " << sigmoid_diff  << std::endl;

			++count;
			score += sigmoid_diff;
		}
	}

	if ( count != 0 ) {
		score = ( score / count ) * n_aa;
	}

	//score *= 3;
	//std::cout << "score(" << pos1 << "," << pos2 << ") = " << score
	//	<< std::endl << std::endl;

	return score;
} // score

} // sequence
} // core
