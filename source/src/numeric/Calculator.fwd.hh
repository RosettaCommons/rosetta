// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking/Calculator.hh
/// @brief  a string-input calculator
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_numeric_Calculator_fwd_hh
#define INCLUDED_numeric_Calculator_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace numeric {

class Calculator;

typedef utility::pointer::shared_ptr<Calculator> CalculatorOP;

} // numeric

#endif // INCLUDED_numeric_Calculator_FWD_HH


