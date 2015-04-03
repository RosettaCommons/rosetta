// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProfSimScoringScheme.hh
/// @brief class definition for a given scoring scheme for an alignment.
/// @details Simply based on comparing single profiles from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ProfSimScoringScheme.hh>

#include <utility/exit.hh>

#include <core/chemical/AA.hh>

#include <string>

#include <utility/vector1.hh>
#include <complex>
#include <map>


namespace core {
namespace sequence {

void ProfSimScoringScheme::initialize_parameters() {
	// this is stupid, refactor this into a database file someday!
	// log-probabilities (base e) for each amino acid calculated from
	// frequencies measured in the SwissProt database on 10/10/08.
	std::map< char, Real > log_p_aa_;
	log_p_aa_['A'] = -2.50775275017596;
	log_p_aa_['C'] = -4.25443370709812;
	log_p_aa_['D'] = -2.91635539230326;
	log_p_aa_['E'] = -2.69728310210688;
	log_p_aa_['F'] = -3.24970340180708;
	log_p_aa_['G'] = -2.6536761536201;
	log_p_aa_['H'] = -3.77884402329716;
	log_p_aa_['I'] = -2.82635555265843;
	log_p_aa_['K'] = -2.83305056323539;
	log_p_aa_['L'] = -2.33558638262581;
	log_p_aa_['M'] = -3.72639555635987;
	log_p_aa_['N'] = -3.20234385542140;
	log_p_aa_['P'] = -3.0416188351623;
	log_p_aa_['Q'] = -3.22656571396179;
	log_p_aa_['R'] = -2.89952136047537;
	log_p_aa_['S'] = -2.70771691574528;
	log_p_aa_['T'] = -2.9259691035378;
	log_p_aa_['V'] = -2.68564236933079;
	log_p_aa_['W'] = -4.51369283589433;
	log_p_aa_['Y'] = -3.53006197206897;

	// initialize log_p_aa values into prior_probs based on the ordering
	// in AA enum
	prior_probs_.resize( log_p_aa_.size() );
	using std::map;
	for ( map< char, Real >::const_iterator it = log_p_aa_.begin(),
			end = log_p_aa_.end();
			it != end; ++it
	) {
			core::chemical::AA aa
				= core::chemical::aa_from_oneletter_code(it->first );
			prior_probs_[ aa ] = std::exp( it->second );
	}
}

Real ProfSimScoringScheme::score(
	SequenceOP seq1,
	SequenceOP seq2,
	Size pos1,
	Size pos2
) {
	SequenceProfileOP prof1
		= SequenceProfileOP( utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq1 ) );
	SequenceProfileOP prof2
		= SequenceProfileOP( utility::pointer::static_pointer_cast< core::sequence::SequenceProfile > ( seq2 ) );

	runtime_assert( pos1 <= prof1->length() );
	runtime_assert( pos2 <= prof2->length() );
	//runtime_assert( (*prof1)[pos1].size() == (*prof2)[pos2].size() );
	runtime_assert( prof1->prof_row( pos1 ).size() == prof2->prof_row(pos2).size() );

	// initialize prior probabilities
	Real divergence_score( 0.0 );
	Real similarity_score( 0.0 );
	Real const base( 2.0 ); // calculate logarithms in base 2

	Size n_aa( prof1->prof_row(pos1).size() );
	//Size n_aa( prof1->alphabet().size() );
	for ( Size i = 1; i <= n_aa; ++i ) {
		// divergence_score is divergence between prof1 and prof2
		Real div_avg( ( prof1->prof_row(pos1)[i] + prof2->prof_row(pos2)[i] ) / 2 );
		divergence_score += 0.5 * prof1->prof_row(pos1)[i]
			* log( prof1->prof_row(pos1)[i] / div_avg ) / log( base );
		divergence_score += 0.5 * prof2->prof_row(pos2)[i]
			* log( prof2->prof_row(pos2)[i] / div_avg ) / log( base );

		// similarity_score is divergence between average of prof1 and prof2 and the prior
		Real prior_prob = prior_probs_[i];
		Real sim_avg( (div_avg + prior_prob) / 2 );
		similarity_score += 0.5 * sim_avg    * log( sim_avg    / prior_prob ) / log( base );
		similarity_score += 0.5 * prior_prob * log( prior_prob / sim_avg    ) / log( base );

		//std::cout << "comparing " << prof1->prof_row(pos1)[i] << " with " << prof2->prof_row(pos2)[i] << std::endl;
	}

	Real score = 0.5 * ( 1 - divergence_score ) * ( 1 + similarity_score );
	//std::cout << "divergence = " << divergence_score
	//	<< ", similarity = " << similarity_score
	//	<< ", score = " << score
	//	<< std::endl;
	return score;
} // score

} // sequence
} // core
