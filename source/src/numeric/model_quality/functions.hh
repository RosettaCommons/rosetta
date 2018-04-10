// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/model_quality/rms.hh
/// @brief
/// @author Jared Adolf-Bryfogle



#ifndef INCLUDED_numeric_model_quality_functions_HH
#define INCLUDED_numeric_model_quality_functions_HH

#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace model_quality {

///@brief Calculate the dihedral angle distance from directional statistics.
///
///@details
/// Metric described in:
///North B, Lehmann A, Dunbrack RL. A new clustering of antibody CDR loop conformations. J Mol Biol 2011; 406:228-256.
numeric::Real
calculate_dihedral_distance( numeric::Real dih1, numeric::Real dih2);


} // end namespace model_quality
} // end namespace numeric

#endif
