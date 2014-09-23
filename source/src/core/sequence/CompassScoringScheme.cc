// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CompassScoringScheme.cc
/// @brief method implementations for CompassScoringScheme class.
/// @author James Thompson

#include <core/types.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/CompassScoringScheme.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>

// AUTO-REMOVED #include <core/chemical/AA.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <iostream>
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

	Real CompassScoringScheme::score(
		SequenceOP seq1,
		SequenceOP seq2,
		Size pos1,
		Size pos2
	) {
		SequenceProfileOP prof1 = SequenceProfileOP( utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq1 ) );
		SequenceProfileOP prof2 = SequenceProfileOP( utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq2 ) );

		runtime_assert( pos1 <= prof1->length() );
		runtime_assert( pos2 <= prof2->length() );
		runtime_assert( prof1->prof_row(pos1).size() == prof2->prof_row(pos2).size() );

		// compare the two profiles using Compass metric, which is defined as:
		// S =   c1 * sum( n(1,i) * log( Q(2,i) / p(i) ) )
		//     + c2 * sum( n(2,i) * log( Q(1,i) / p(i) ) )
		// the terms are defined as such:
		// i - column of the profile to be evaluated, represents a residue
		// p(i) - prior probability of residue i in all sequences
		// n(1,i) - the effective number of sequences from profile 1 with residue i
		// Q(2,i) - the estimated frequency of residue i in sequence 2
		// n(2,i) - the effective number of sequences from profile 2 with residue i
		// Q(1,i) - the estimated frequency of residue i in sequence 1
		// c1     - normalization constant for profile 1
		// c2     - normalization constant for profile 2

		// naive. better to calculate using Neff from PSIC weighting.
		//Real const c1(0.5), c2(0.5);
		//c1 = 0.5;
		//c2 = 0.5;

		// calculate normalization constants c1 and c2
		Size n_aa1( prof1->prof_row(pos1).size() );
		Size n_aa2( prof2->prof_row(pos2).size() );
		Real score( 0.0 );
		runtime_assert( n_aa1 == n_aa2 );
		//for ( Size i = 1; i <= n_aa; ++i ) {

		//}

		return score;
	} // score


} // sequence
} // core
