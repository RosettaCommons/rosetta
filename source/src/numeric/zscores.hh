// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/util.hh
/// @brief small bundle of utilities for dealing with zscores.
/// @author Jared Adolf-Bryfogle

#ifndef INCLUDED_numeric_zscores_hh
#define INCLUDED_numeric_zscores_hh

#include <numeric/types.hh>
#include <numeric/numeric.functions.hh>
#include <utility/vector1.hh>
#include <utility/numbers.hh>

#include <numeric/MathNTensorBase.fwd.hh>

#include <limits>
#include <cmath>
#include <algorithm>
#include <map>

namespace numeric {


///@brief Calculate a Z-score from a set of data.
///  Real i_zscore =  (input_v[i]-mean)/stdev;
///
///@author Ray Wang (wangyr@uw.edu)
/// Negating flips the zscore (i_zscore = -1*i_zscore)
void
calc_zscore(
	std::map< Size, Real > const & input_v,
	std::map< Size, Real > & zscore_v,
	bool negating = false
);

} // numeric

#endif
