// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/optimize_weights/Arithmetic.fwd.hh
/// @brief  Parse tree for arithmetic operations class forward declarations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_optimize_weights_Arithmetic_fwd_hh
#define INCLUDED_protocols_optimize_weights_Arithmetic_fwd_hh

/// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace optimize_weights {

///// TOKEN CLASSES /////
class Token;
typedef utility::pointer::owning_ptr< Token > TokenOP;
typedef utility::pointer::owning_ptr< Token const > TokenCOP;

class LiteralToken;
typedef utility::pointer::owning_ptr< LiteralToken > LiteralTokenOP;
typedef utility::pointer::owning_ptr< LiteralToken const > LiteralTokenCOP;

class VariableToken;
typedef utility::pointer::owning_ptr< VariableToken > VariableTokenOP;
typedef utility::pointer::owning_ptr< VariableToken const > VariableTokenCOP;

class FunctionToken;
typedef utility::pointer::owning_ptr< FunctionToken > FunctionTokenOP;
typedef utility::pointer::owning_ptr< FunctionToken const > FunctionTokenCOP;

class SimpleToken;
typedef utility::pointer::owning_ptr< SimpleToken > SimpleTokenOP;
typedef utility::pointer::owning_ptr< SimpleToken const > SimpleTokenCOP;


//// TOKEN CONTAINER /////

class TokenSet;
typedef utility::pointer::owning_ptr< TokenSet > TokenSetOP;
typedef utility::pointer::owning_ptr< TokenSet const > TokenSetCOP;


///// SCANNER CLASS /////

class ArithmeticScanner;
typedef utility::pointer::owning_ptr< ArithmeticScanner > ArithmeticScannerOP;
typedef utility::pointer::owning_ptr< ArithmeticScanner const > ArithmeticScannerCOP;


///// ABSTRACT SYNTAX TREE CLASSES /////

class ArithmeticASTNode;
typedef utility::pointer::owning_ptr< ArithmeticASTNode > ArithmeticASTNodeOP;
typedef utility::pointer::owning_ptr< ArithmeticASTNode const > ArithmeticASTNodeCOP;

class ArithmeticASTExpression;
typedef utility::pointer::owning_ptr< ArithmeticASTExpression > ArithmeticASTExpressionOP;
typedef utility::pointer::owning_ptr< ArithmeticASTExpression const > ArithmeticASTExpressionCOP;

class ArithmeticASTFunction;
typedef utility::pointer::owning_ptr< ArithmeticASTFunction > ArithmeticASTFunctionOP;
typedef utility::pointer::owning_ptr< ArithmeticASTFunction const > ArithmeticASTFunctionCOP;

class ArithmeticASTTerm;
typedef utility::pointer::owning_ptr< ArithmeticASTTerm > ArithmeticASTTermOP;
typedef utility::pointer::owning_ptr< ArithmeticASTTerm const > ArithmeticASTTermCOP;

class ArithmeticASTFactor;
typedef utility::pointer::owning_ptr< ArithmeticASTFactor > ArithmeticASTFactorOP;
typedef utility::pointer::owning_ptr< ArithmeticASTFactor const > ArithmeticASTFactorCOP;

class ArithmeticASTValue;
typedef utility::pointer::owning_ptr< ArithmeticASTValue > ArithmeticASTValueOP;
typedef utility::pointer::owning_ptr< ArithmeticASTValue const > ArithmeticASTValueCOP;

class ArithmeticASTRestTerm;
typedef utility::pointer::owning_ptr< ArithmeticASTRestTerm > ArithmeticASTRestTermOP;
typedef utility::pointer::owning_ptr< ArithmeticASTRestTerm const > ArithmeticASTRestTermCOP;

class ArithmeticASTRestExpression;
typedef utility::pointer::owning_ptr< ArithmeticASTRestExpression > ArithmeticASTRestExpressionOP;
typedef utility::pointer::owning_ptr< ArithmeticASTRestExpression const > ArithmeticASTRestExpressionCOP;


///// AST Visitors /////

class ASTVisitor;
typedef utility::pointer::owning_ptr< ASTVisitor > ASTVisitorOP;
typedef utility::pointer::owning_ptr< ASTVisitor const > ASTVisitorCOP;

class ASTPrinter;
typedef utility::pointer::owning_ptr< ASTPrinter > ASTPrinterOP;
typedef utility::pointer::owning_ptr< ASTPrinter const > ASTPrinterCOP;


//// EXPRESSION CLASSES ////

class Expression;
typedef utility::pointer::owning_ptr< Expression > ExpressionOP;
typedef utility::pointer::owning_ptr< Expression const > ExpressionCOP;

class LiteralExpression;
typedef utility::pointer::owning_ptr< LiteralExpression > LiteralExpressionOP;
typedef utility::pointer::owning_ptr< LiteralExpression const > LiteralExpressionCOP;

class VariableExpression;
typedef utility::pointer::owning_ptr< VariableExpression > VariableExpressionOP;
typedef utility::pointer::owning_ptr< VariableExpression const > VariableExpressionCOP;

class UnaryExpression;
typedef utility::pointer::owning_ptr< UnaryExpression > UnaryExpressionOP;
typedef utility::pointer::owning_ptr< UnaryExpression const > UnaryExpressionCOP;

class BinaryExpression;
typedef utility::pointer::owning_ptr< BinaryExpression > BinaryExpressionOP;
typedef utility::pointer::owning_ptr< BinaryExpression const > BinaryExpressionCOP;

class SquarerootExpression;
typedef utility::pointer::owning_ptr< SquarerootExpression > SquarerootExpressionOP;
typedef utility::pointer::owning_ptr< SquarerootExpression const > SquarerootExpressionCOP;

class AddExpression;
typedef utility::pointer::owning_ptr< AddExpression > AddExpressionOP;
typedef utility::pointer::owning_ptr< AddExpression const > AddExpressionCOP;

class SubtractExpression;
typedef utility::pointer::owning_ptr< SubtractExpression > SubtractExpressionOP;
typedef utility::pointer::owning_ptr< SubtractExpression const > SubtractExpressionCOP;

class MultiplyExpression;
typedef utility::pointer::owning_ptr< MultiplyExpression > MultiplyExpressionOP;
typedef utility::pointer::owning_ptr< MultiplyExpression const > MultiplyExpressionCOP;

class DivideExpression;
typedef utility::pointer::owning_ptr< DivideExpression > DivideExpressionOP;
typedef utility::pointer::owning_ptr< DivideExpression const > DivideExpressionCOP;

class MaxExpression;
typedef utility::pointer::owning_ptr< MaxExpression > MaxExpressionOP;
typedef utility::pointer::owning_ptr< MaxExpression const > MaxExpressionCOP;

class MinExpression;
typedef utility::pointer::owning_ptr< MinExpression > MinExpressionOP;
typedef utility::pointer::owning_ptr< MinExpression const > MinExpressionCOP;

class MetaMaxExpression;
typedef utility::pointer::owning_ptr< MetaMaxExpression > MetaMaxExpressionOP;
typedef utility::pointer::owning_ptr< MetaMaxExpression const > MetaMaxExpressionCOP;

class MetaMinExpression;
typedef utility::pointer::owning_ptr< MetaMinExpression > MetaMinExpressionOP;
typedef utility::pointer::owning_ptr< MetaMinExpression const > MetaMinExpressionCOP;

}
}

#endif

