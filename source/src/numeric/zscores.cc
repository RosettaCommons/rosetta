// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/util.cc
/// @brief small bundle of utilities for dealing with zscores.
/// @author Jared Adolf-Bryfogle

#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <numeric/MathNTensorBase.hh>
#include <numeric/MathNTensor.hh>

#include <algorithm>
#include <map>
#include <cmath>

namespace numeric {


void
calc_zscore(
	std::map< Size, Real > const & input_v,
	std::map< Size, Real > & zscore_v,
	bool negating
){
	
	Real sum=0, sq_sum=0;
	Size nres = input_v.size();

	for ( auto score_pair: input_v ) {
		sum    += score_pair.second;
		sq_sum += score_pair.second * score_pair.second;
	}
	Real mean  = sum/nres;
	Real stdev = std::sqrt( sq_sum/nres - mean * mean );
	
	//std::cout << mean << " " << stdev << std::endl;
	
	for ( auto score_pair : input_v ) {
		
		//std::cout << "rscore: " << score_pair.second << std::endl;
		Real i_zscore =  (score_pair.second - mean)/stdev;
		
		//std::cout << "iscore: " << i_zscore << std::endl;
		if ( negating ) i_zscore = -1*i_zscore;

		zscore_v[ score_pair.first] = i_zscore;
	}
}

} // numeric
