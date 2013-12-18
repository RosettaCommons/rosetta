// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
typedef utility::pointer::owning_ptr< VectorExpression > VectorExpressionOP;
typedef utility::pointer::owning_ptr< VectorExpression const > VectorExpressionCOP;

class VariableVectorExpression;
typedef utility::pointer::owning_ptr< VariableVectorExpression > VariableVectorExpressionOP;
typedef utility::pointer::owning_ptr< VariableVectorExpression const > VariableVectorExpressionCOP;

class IterativeVectorExpression;
typedef utility::pointer::owning_ptr< IterativeVectorExpression > IterativeVectorExpressionOP;
typedef utility::pointer::owning_ptr< IterativeVectorExpression const > IterativeVectorExpressionCOP;

class VectorFunction;
typedef utility::pointer::owning_ptr< VectorFunction > VectorFunctionOP;
typedef utility::pointer::owning_ptr< VectorFunction const > VectorFunctionCOP;

class VMax;
typedef utility::pointer::owning_ptr< VMax > VMaxOP;
typedef utility::pointer::owning_ptr< VMax const > VMaxCOP;

class VMin;
typedef utility::pointer::owning_ptr< VMin > VMinOP;
typedef utility::pointer::owning_ptr< VMin const > VMinCOP;

class PowExpression;
typedef utility::pointer::owning_ptr< PowExpression > PowExpressionOP;
typedef utility::pointer::owning_ptr< PowExpression const > PowExpressionCOP;

class ExpExpression;
typedef utility::pointer::owning_ptr< ExpExpression > ExpExpressionOP;
typedef utility::pointer::owning_ptr< ExpExpression const > ExpExpressionCOP;

class LnExpression;
typedef utility::pointer::owning_ptr< LnExpression > LnExpressionOP;
typedef utility::pointer::owning_ptr< LnExpression const > LnExpressionCOP;

class InSetExpression;
typedef utility::pointer::owning_ptr< InSetExpression > InSetExpressionOP;
typedef utility::pointer::owning_ptr< InSetExpression const > InSetExpressionCOP;

class VectorExpressionCreator;
typedef utility::pointer::owning_ptr< VectorExpressionCreator > VectorExpressionCreatorOP;
typedef utility::pointer::owning_ptr< VectorExpressionCreator const > VectorExpressionCreatorCOP;

class SurrogateVariableExpression;
typedef utility::pointer::owning_ptr< SurrogateVariableExpression > SurrogateVariableExpressionOP;
typedef utility::pointer::owning_ptr< SurrogateVariableExpression const > SurrogateVariableExpressionCOP;

class FileContentsMap;
typedef utility::pointer::owning_ptr< FileContentsMap > FileContentsMapOP;
typedef utility::pointer::owning_ptr< FileContentsMap const > FileContentsMapCOP;

class DynamicAggregateFunction;

typedef utility::pointer::owning_ptr< DynamicAggregateFunction > DynamicAggregateFunctionOP;
typedef utility::pointer::owning_ptr< DynamicAggregateFunction const > DynamicAggregateFunctionCOP;

class DynamicAggregateFunctionDriver;

typedef utility::pointer::owning_ptr< DynamicAggregateFunctionDriver > DynamicAggregateFunctionDriverOP;
typedef utility::pointer::owning_ptr< DynamicAggregateFunctionDriver const > DynamicAggregateFunctionDriverCOP;

class EntityFuncExpressionCreator;

class EntityFunc;
typedef utility::pointer::owning_ptr< EntityFunc > EntityFuncOP;
typedef utility::pointer::owning_ptr< EntityFunc const > EntityFuncCOP;

}
}

#endif
