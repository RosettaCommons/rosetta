// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/OptEMultifunc.fwd.hh
/// @brief  core::optimization::OptEMultifunc forward declarations
/// @author Jim Havranek


#ifndef INCLUDED_protocols_optimize_weights_OptEMultifunc_fwd_hh
#define INCLUDED_protocols_optimize_weights_OptEMultifunc_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace optimize_weights {

class OptEMultifunc;
typedef utility::pointer::shared_ptr< OptEMultifunc > OptEMultifuncOP;
typedef utility::pointer::shared_ptr< OptEMultifunc const > OptEMultifuncCOP;
typedef utility::pointer::weak_ptr< OptEMultifunc > OptEMultifuncAP;
typedef utility::pointer::weak_ptr< OptEMultifunc const > OptEMultifuncCAP;

class WrapperOptEMultifunc;
typedef utility::pointer::shared_ptr< WrapperOptEMultifunc > WrapperOptEMultifuncOP;
typedef utility::pointer::shared_ptr< WrapperOptEMultifunc const > WrapperOptEMultifuncCOP;
typedef utility::pointer::weak_ptr< WrapperOptEMultifunc > WrapperOptEMultifuncAP;
typedef utility::pointer::weak_ptr< WrapperOptEMultifunc const > WrapperOptEMultifuncCAP;

class OptEVariableExpression;
typedef utility::pointer::shared_ptr< OptEVariableExpression > OptEVariableExpressionOP;
typedef utility::pointer::shared_ptr< OptEVariableExpression const > OptEVariableExpressionCOP;


} // namespace optimize_weights
} // namespace protocols


#endif // INCLUDED_protocols_optimize_weights_OptEMultifunc_FWD_HH
