// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random.functions.cc
/// @brief  Random number functions
/// @author Original author unknown
///
/// @note These are the rosetta++ functions that will be replaced by the numeric::random package to
/// provide decoupled random generator objects for solution stability of test cases
///
/// @note This is purposely NOT in the numeric namespace to simplify migration


// Unit headers
#include <numeric/random.functions.hh>

#include <numeric/random/random.hh>

// Project headers
// #include <basic/options/option.hh>

// Utility Headers
#include <utility/exit.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/Time_Date.hh>

// C++ Headers
#include <cmath>
#include <ctime>
#include <iostream>


