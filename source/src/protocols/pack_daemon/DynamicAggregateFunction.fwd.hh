// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/MultistateAggregateFunction.fwd.hh
/// @brief  forward declaration of class MultistateAggregateFunction
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_DynamicAggregateFunction_fwd_hh
#define INCLUDED_protocols_pack_daemon_DynamicAggregateFunction_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace pack_daemon {

class VectorExpression;
typedef utility::pointer::shared_ptr< VectorExpression > VectorExpressionOP;
typedef utility::pointer::shared_ptr< VectorExpression const > VectorExpressionCOP;

class VariableVectorExpression;
typedef utility::pointer::shared_ptr< VariableVectorExpression > VariableVectorExpressionOP;
typedef utility::pointer::shared_ptr< VariableVectorExpression const > VariableVectorExpressionCOP;

class IterativeVectorExpression;
typedef utility::pointer::shared_ptr< IterativeVectorExpression > IterativeVectorExpressionOP;
typedef utility::pointer::shared_ptr< IterativeVectorExpression const > IterativeVectorExpressionCOP;

class VectorFunction;
typedef utility::pointer::shared_ptr< VectorFunction > VectorFunctionOP;
typedef utility::pointer::shared_ptr< VectorFunction const > VectorFunctionCOP;

class VMax;
typedef utility::pointer::shared_ptr< VMax > VMaxOP;
typedef utility::pointer::shared_ptr< VMax const > VMaxCOP;

class VMin;
typedef utility::pointer::shared_ptr< VMin > VMinOP;
typedef utility::pointer::shared_ptr< VMin const > VMinCOP;

class PowExpression;
typedef utility::pointer::shared_ptr< PowExpression > PowExpressionOP;
typedef utility::pointer::shared_ptr< PowExpression const > PowExpressionCOP;

class ExpExpression;
typedef utility::pointer::shared_ptr< ExpExpression > ExpExpressionOP;
typedef utility::pointer::shared_ptr< ExpExpression const > ExpExpressionCOP;

class LnExpression;
typedef utility::pointer::shared_ptr< LnExpression > LnExpressionOP;
typedef utility::pointer::shared_ptr< LnExpression const > LnExpressionCOP;

class InSetExpression;
typedef utility::pointer::shared_ptr< InSetExpression > InSetExpressionOP;
typedef utility::pointer::shared_ptr< InSetExpression const > InSetExpressionCOP;

class VectorExpressionCreator;
typedef utility::pointer::shared_ptr< VectorExpressionCreator > VectorExpressionCreatorOP;
typedef utility::pointer::shared_ptr< VectorExpressionCreator const > VectorExpressionCreatorCOP;

class SurrogateVariableExpression;
typedef utility::pointer::shared_ptr< SurrogateVariableExpression > SurrogateVariableExpressionOP;
typedef utility::pointer::shared_ptr< SurrogateVariableExpression const > SurrogateVariableExpressionCOP;

class FileContentsMap;
typedef utility::pointer::shared_ptr< FileContentsMap > FileContentsMapOP;
typedef utility::pointer::shared_ptr< FileContentsMap const > FileContentsMapCOP;

class DynamicAggregateFunction;

typedef utility::pointer::shared_ptr< DynamicAggregateFunction > DynamicAggregateFunctionOP;
typedef utility::pointer::shared_ptr< DynamicAggregateFunction const > DynamicAggregateFunctionCOP;

class DynamicAggregateFunctionDriver;

typedef utility::pointer::shared_ptr< DynamicAggregateFunctionDriver > DynamicAggregateFunctionDriverOP;
typedef utility::pointer::shared_ptr< DynamicAggregateFunctionDriver const > DynamicAggregateFunctionDriverCOP;

class EntityFuncExpressionCreator;

class EntityFunc;
typedef utility::pointer::shared_ptr< EntityFunc > EntityFuncOP;
typedef utility::pointer::shared_ptr< EntityFunc const > EntityFuncCOP;

}
}

#endif
