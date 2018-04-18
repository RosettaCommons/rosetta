// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/model_quality/functions.hh
/// @brief
/// @author Jared Adolf-Bryfogle


// Rosetta Headers

// numeric libraries
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/types.hh>

#include <utility/vector1.hh>


namespace numeric {
namespace model_quality {

numeric::Real
calculate_dihedral_distance( numeric::Real dih1, numeric::Real dih2 ){
	return ( 2 * ( 1- cos ( ( dih1- dih2 ) * numeric::NumericTraits<numeric::Real>::pi() /180) ) );
}


} // rms
} // numeric
