// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/PackingScore.cc
///
/// @brief
/// @author will sheffler


// Unit header or inline function header
#include <core/scoring/packstat/PackingScore.hh>

#include <cmath>
#include <iostream>


namespace core {
namespace scoring {
namespace packstat {

using core::Real;
using utility::vector1;

// bool print = false;

std::ostream & operator<< ( std::ostream & out, PackingScoreResData const & dat ) {
	for ( Size i =1; i <= dat.nrad(); ++i ) {
		for ( Size j =1; j <= dat.npr(); ++j ) {
			out << dat.msa(i,j) << " ";
		}
		out << "     ";
	}
	return out;
}


Real PackingScore::score( PackingScoreResDataCOP dat ) const {
	debug_assert( dat->npr() == npr() && dat->nrad() == nrad() );
	Real score = 0;
	for ( Size i =1; i <= nrad(); ++i ) {
		for ( Size j =1; j <= npr(); ++j ) {
			// std::cerr << "DAT " << i << " " << j << " " << dat->msa(i,j) << std::endl;
			score += ( dat->msa(i,j) - center(i,j) ) * weight(i,j);
		}
	}
	score -= rho();
	// if(print) std::cerr << " DV " << score;
	if ( compprob() ) {
		score = 1.0 - (1.0 / (1.0 + exp( probA_ * score + probB_ ) ));
	} else {
		score = probA_ * score + probB_;
	}
	// if(print) std::cerr << " SC " << score << std::endl;
	return score;
}


Real PackingScore::score( utility::vector1<PackingScoreResDataCOP> dats ) const {
	Real sc = 0.0/*, psc = 0.0*/;
	for ( Size di = 1; di <= dats.size(); ++di ) {
		for ( Size i =1; i <= nrad(); ++i ) {
			for ( Size j =1; j <= npr(); ++j ) {
				// std::cerr << i << " " << j << " DAT " << dats[di]->msa(i,j) << std::endl;
				sc += ( dats[di]->msa(i,j) - center(i,j) ) * weight(i,j);
			}
		}
		// if( di == 1 ) {
		//  std::cerr << "di " << di;
		//  print = true;
		// } else {
		//  print = false;
		// }
		//
		// psc += score( dats[di] );
		// std::cerr << "score: " << di << " " << score(dats[di]) << std::endl;
		// std::exit(-1);
	}
	sc = sc / dats.size();
	sc = sc - rho();
	// std::cout << "SCORE " << sc << std::endl;
	if ( compprob() ) {
		sc = 1.0 - (1.0 / (1.0 + exp( probA_ * sc + probB_ ) ) );
	} else {
		sc = probA_ * sc + probB_;
	}
	// std::cerr << "avg: " << psc / (float)dats.size() << " " << sc << std::endl;
	return sc;
}


} // namespace packstat
} // namespace scoring
} // namespace core
