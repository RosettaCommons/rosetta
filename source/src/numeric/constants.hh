// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/constants.hh
/// @brief  Common numeric constants in varying precisions
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @note   The 'constants' namespace and the namespaces within it
//          do not correspond to any package


#ifndef INCLUDED_numeric_constants_hh
#define INCLUDED_numeric_constants_hh


#include <numeric/types.hh>


namespace numeric {
namespace constants {


// float
namespace f {

typedef  float  Type;

extern Type const zero;
extern Type const one;
extern Type const two;
extern Type const pi;
extern Type const pi_2;              // 2 * pi
extern Type const pi_over_2;         // pi / 2
extern Type const pi_over_3;         // pi / 3
extern Type const pi_2_over_3;       // ( 2 * pi ) / 3
extern Type const pi_over_180;       // pi / 180
extern Type const degrees_to_radians;
extern Type const deg2rad;
extern Type const radians_to_degrees;
extern Type const rad2deg;

} // namespace f


// double
namespace d {

typedef  double  Type;

extern Type const zero;
extern Type const one;
extern Type const two;
extern Type const pi;
extern Type const pi_2;             // 2 * pi
extern Type const pi_over_2;        // pi / 2
extern Type const pi_over_3;        // pi / 3
extern Type const pi_2_over_3;      // ( 2 * pi ) / 3
extern Type const pi_over_180;      // pi / 180
extern Type const degrees_to_radians;
extern Type const deg2rad;
extern Type const radians_to_degrees;
extern Type const rad2deg;

} // namespace d


// long double
namespace ld {

typedef  long double  Type;

extern Type const zero;
extern Type const one;
extern Type const two;
extern Type const pi;
extern Type const pi_2;             // 2 * pi
extern Type const pi_over_2;        // pi / 2
extern Type const pi_over_3;        // pi / 3
extern Type const pi_2_over_3;      // ( 2 * pi ) / 3
extern Type const pi_over_180;      // pi / 180
extern Type const degrees_to_radians;
extern Type const deg2rad;
extern Type const radians_to_degrees;
extern Type const rad2deg;

} // namespace ld


// Real
namespace r {

typedef  Real  Type;

extern Type const zero;
extern Type const one;
extern Type const two;
extern Type const pi;
extern Type const pi_2;             // 2 * pi
extern Type const pi_over_2;        // pi / 2
extern Type const pi_over_3;        // pi / 3
extern Type const pi_2_over_3;      // ( 2 * pi ) / 3
extern Type const pi_over_180;      // pi / 180
extern Type const degrees_to_radians;
extern Type const deg2rad;
extern Type const radians_to_degrees;
extern Type const rad2deg;

} // namespace r


} // namespace constants
} // namespace numeric


#endif // INCLUDED_numeric_constants_HH
