// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>

#include <core/sequence/Sequence.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <ObjexxFCL/format.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

///////////////////////////////////////////////////////////////////////////////
using std::map;
using std::min;
using std::max;
using std::exp;
using std::pair;
using std::make_pair;
using core::Real;
using core::Size;
using std::string;
using utility::vector1;

using namespace core::sequence;

MultipleSequenceAlignment random_alignment_move( MultipleSequenceAlignment const & msa ) {
	MultipleSequenceAlignment new_msa = msa;
	// pick a random sequence
	Size sequence_idx = numeric::random::random_range( 1, (int) new_msa.size() );

	// randomly pick a position in the alignment.
	Size position = numeric::random::random_range( 5, (int) new_msa.sequence(sequence_idx).size() - 5 );

	// if there's a gap here, either delete it or extend it
	// Real delete_prob = 0.5;
	// if ( new_msa.sequence(sequence_idx).is_gap(position) ) {
	// 	if ( numeric::random::uniform() > delete_prob ) {
	//
	// 	}
	// }
	// std::cout << "inserting gap into sequence " << sequence_idx
	//  					<< ", position " << position << std::endl;

	new_msa.insert_gap_into_sequence( sequence_idx, position );
	return new_msa;
} // random_alignment_move

Real blosum_score_alignment( MultipleSequenceAlignment & msa ) {
	static bool init( false );
	using namespace core::chemical;
	// static map< string, map< string, Real > > scoring_matrix;

	vector1< Real > dummy( num_aa_types, 0.0 );
	vector1< vector1< Real > > scoring_matrix( num_aa_types, dummy ); // indexed by entries in AA enum

	if ( !init ) {
		utility::io::izstream input( "BLOSUM62" );
		if ( !input ) {
			utility_exit_with_message( "ERROR: Unable to open BLOSUM62!" );
		}

		string line;
		vector1< AA > order;
		while( getline( input, line ) ) {
			if ( line.substr(0,1) == "#" ) continue; // skip comments

			std::istringstream line_stream( line );
			if ( line.substr(0,1) == " " ) { // header line
				while ( !line_stream.fail() ) {
					char aa;
					line_stream >> aa;
					if ( oneletter_code_specifies_aa( aa ) )
						order.push_back( aa_from_oneletter_code(aa) );
				}
			} else {
				char aa_name;
				line_stream >> aa_name;
				if ( !oneletter_code_specifies_aa( aa_name ) ) continue; // skip non-AA's
				AA aa = aa_from_oneletter_code(aa_name);

				for ( vector1< AA >::const_iterator it = order.begin(), end = order.end();
							it != end; ++it
				) {
					Real score;
					line_stream >> score;
					scoring_matrix[ aa ][ *it ] = score;

					if ( line_stream.fail() ) {
						string message = "Error reading line" + line + '\n';
						utility_exit_with_message( message );
					}
				}
			}
		} // while( getline( input, line ) )

		init = true;
	} // if ( !init )

	Real score( 0 );
	Real gap_penalty = -6;

	for ( Size p = 1; p <= msa.sequence(1).size(); ++p ) {
		for ( Size i = 1; i <= msa.size(); ++i ) {
			for ( Size j = i+1; j <= msa.size(); ++j ) {
				char char1 = msa.sequence(i)[p];
				char char2 = msa.sequence(j)[p];
				if ( msa.sequence(i).is_gap( p ) || msa.sequence(j).is_gap( p ) ) {
					score += gap_penalty;
				} else {
					Real local_score
						= scoring_matrix[ aa_from_oneletter_code(char1) ][ aa_from_oneletter_code(char2) ];
					// if ( char1 == char2 ) local_score = 4;
					// else local_score = 1;
					score += local_score;
					// std::cout << "score(" << char1 << "," << char2 << ") = " << local_score << std::endl;
				}
			} // j
		} // i
	} // p


	score = -1 * score;
	msa.score( score );
	std::cout << msa << std::endl;

	return score;
} // blosum_score_alignment


Real score_alignment( MultipleSequenceAlignment & msa ) {
	return blosum_score_alignment( msa );
} // score_alignment

Real get_temp(
	Size current_iter,
	Size /*max_iter*/
) {
	// Size offset = 1; // to keep temperature from going absolutely to zero.
	// return max_iter - current_iter + offset;
	// stupid exponential cooling scheme!
	// Real alpha = 1 - ( 1 / max_iter );
	Real alpha = 0.95;
	Real temp = 1000;
	for ( Size i = 1; i <= current_iter; ++i ) {
		temp = temp * alpha;
	}
	return temp;
}

int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );

	Sequence seq1(
		// "VDKF",
		// "VDDNKFAGLMLVGTIEILHDRASKEMLWTDGCEIYYPLGIDDPDYTALCFTAEWGNYYRH",
		"IDEKFLIESNELVESSKIVMVGTNGENGYPNIKAMMRLKHDGLKKFWLSTNTSTRMVERLKKNNKICLYFVDDNKFAGLMLVGTIEILHDRASKEMLWTDGCEIYYPLGIDDPDYTALCFTAEWGNYYRHLKNITFKIDEIY",
		"t380_"
	);

	Sequence seq2(
		// "VDDQQQNKF",
		// "QE--KGDSVALXGEVEVVTDEKLKQELWQDWFIEHFPGGPTDPGYVLLKFTANHATYWIE",
		"TKTMKEKAVELLQKCEVVTLASVNKEGYPRPVPMSKIAAEGISTIWMSTGADSLKTIDFLSNPKAGLCFQEKGDSVALMGEVEVVTDEKLKQELWQDWFIEHFPGGPTDPGYVLLKFTANHATYWIEGTFIHKKL",
		"2fhqA"
	);

	vector1< MultipleSequenceAlignment > results;

	for ( Size run_count = 1; run_count <= 100; ++run_count ) {
		MultipleSequenceAlignment msa, last_accepted_msa, best_msa;

		msa.add_sequence( seq1 );
		msa.add_sequence( seq2 );

		score_alignment( msa );

		std::cout << msa << std::endl;

		best_msa = msa;
		last_accepted_msa = msa;

		Size nmoves    = 2;
		Size n_iter    = 10000;
		Size n_accepts = 0;
		for ( Size i = 1; i <= n_iter; ++i ) {
			// for ( Size j = 1; j <= 80; ++j ) {
			// 	std::cout << '*';
			// }
			// std::cout << std::endl;

			// std::cout << "iteration: " << i << std::endl;
			for ( Size k = 1; k <= nmoves; ++k )
				msa = random_alignment_move( msa );
			score_alignment( msa );

			std::cout << "scores(current,last_accepted,best) = ("
								<< msa.score() << ","
								<< last_accepted_msa.score() << ","
								<< best_msa.score() << ")"
								<< std::endl;
			std::cout << msa << std::endl;

			if ( msa.score() < best_msa.score() ) {
				std::cout << "accept (best)" << std::endl;
				std::cout << msa << std::endl;
				best_msa          = msa;
				last_accepted_msa = best_msa;
				++n_accepts;

				std::cout << "best_msa: " << best_msa << std::endl;
			} else {
				// score is worse (higher) than last_accepted_score.
				Real const temperature = get_temp( n_accepts, n_iter );
				Real const boltz_factor = ( last_accepted_msa.score() - msa.score() ) / temperature;
				Real const dice_roll    = numeric::random::uniform();
				Real const probability  = exp( min( 40.0, max( -40.0, boltz_factor ) ) );

				std::cout << "trying (" << last_accepted_msa.score() << " -> " << msa.score()
									<< ") with probability " << probability << std::endl;

				if ( dice_roll >= probability ) {
					// std::cout << "rejected" << std::endl;
					msa = last_accepted_msa;
				} else {
					// std::cout << "accept (lucked out with probability "
					// 					<< dice_roll << " <= " << probability	<<  ")"
					// 					<< " temp = " << temperature
					// 					<< std::endl;
					last_accepted_msa = msa;
					++n_accepts;
				}
			} // else

			// std::cout << msa << std::endl;
		} // for n_iter

		std::cout << "done (" << n_accepts << " accepts)" << std::endl;
		results.push_back( best_msa );
	} // k


	for ( vector1< MultipleSequenceAlignment >::iterator
				it = results.begin(), end = results.end();
	 			it != end; ++it
	) {
		for ( Size j = 1; j <= 80; ++j )
			std::cout << '*';
		std::cout << std::endl;

		std::cout << *it << std::endl;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
