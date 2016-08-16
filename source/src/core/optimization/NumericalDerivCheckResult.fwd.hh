// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/NumericalDerivCheckResult.fwd.hh
/// @brief  Forward declaration for nuerical derivative check results classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_optimization_NumericalDerivCheckResult_fwd_hh
#define INCLUDED_core_optimization_NumericalDerivCheckResult_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace optimization {

class DOF_DataPoint;
class SimpleDerivCheckResult;
class NumDerivCheckData;
class NumericalDerivCheckResult;


typedef utility::pointer::shared_ptr< SimpleDerivCheckResult > SimpleDerivCheckResultOP;
typedef utility::pointer::shared_ptr< SimpleDerivCheckResult const > SimpleDerivCheckResultCOP;

typedef utility::pointer::shared_ptr< NumDerivCheckData > NumDerivCheckDataOP;
typedef utility::pointer::shared_ptr< NumDerivCheckData const > NumDerivCheckDataCOP;

typedef utility::pointer::shared_ptr< NumericalDerivCheckResult > NumericalDerivCheckResultOP;
typedef utility::pointer::shared_ptr< NumericalDerivCheckResult const > NumericalDerivCheckResultCOP;

} // namespace optimization
} // namespace core


#endif
