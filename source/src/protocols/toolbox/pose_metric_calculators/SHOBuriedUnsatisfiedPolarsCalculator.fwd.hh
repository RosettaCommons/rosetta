// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief Definition of types SHOBuriedUnsatisfiedPolarsCalculatorOP and
/// 	SHOBuriedUnsatisfiedPolarsCalculatorCOP
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)


#ifndef INCLUDED_FWD_SHOBuriedUnsatisfiedPolarsCalculator_hh
#define INCLUDED_FWD_SHOBuriedUnsatisfiedPolarsCalculator_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class SHOBuriedUnsatisfiedPolarsCalculator;
typedef utility::pointer::shared_ptr< SHOBuriedUnsatisfiedPolarsCalculator > SHOBuriedUnsatisfiedPolarsCalculatorOP;
typedef utility::pointer::shared_ptr< SHOBuriedUnsatisfiedPolarsCalculator const > SHOBuriedUnsatisfiedPolarsCalculatorCOP;

} // pose_metrics_calculatorss
} // toolbox
} // protocols

#endif
