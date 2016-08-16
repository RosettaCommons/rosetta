// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/buns/BuriedUnsatisfiedPolarsCalculator2.fwd.hh
/// @brief
/// @author Kevin Houlihan


#ifndef INCLUDED_devel_buns_BuriedUnsatisfiedPolarsCalculator2_fwd_hh
#define INCLUDED_devel_buns_BuriedUnsatisfiedPolarsCalculator2_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace buns {

class BuriedUnsatisfiedPolarsCalculator2;
enum CALCULATORMODE { energy, geometry, check_buns1 };

typedef utility::pointer::shared_ptr< BuriedUnsatisfiedPolarsCalculator2> BuriedUnsatisfiedPolarsCalculator2OP;
typedef utility::pointer::shared_ptr< BuriedUnsatisfiedPolarsCalculator2 const > BuriedUnsatisfiedPolarsCalculator2COP;

} // namespace toolbox
} // namespace protocols

#endif
