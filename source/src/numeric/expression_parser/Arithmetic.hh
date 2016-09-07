// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/expression_parser/Arithmetic.hh
/// @brief  Parse tree for arithmetic operations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_numeric_expression_parser_Arithmetic_HH
#define INCLUDED_numeric_expression_parser_Arithmetic_HH

/// Unit headers
#include <numeric/expression_parser/Arithmetic.fwd.hh>

/// Numeric headers
#include <numeric/types.hh>

/// Utility headers
#include <utility/pointer/ReferenceCount.hh>

/// C++ headers
#include <sstream>
#include <list>
#include <map>

#include <utility/vector1.hh>


namespace numeric {
namespace expression_parser {

enum TokenType {
	INVALID_TOKEN_TYPE,
	LITERAL,
	VARIABLE,
	FUNCTION,
	COMMA,
	LEFT_PAREN,
	RIGHT_PAREN,
	PLUS_SYMBOL,
	SUBTRACT_SYMBOL,
	MULTIPLY_SYMBOL,
	DIVIDE_SYMBOL//,
	//TERMINAL_SYMBOL
};

std::string
token_type_name( TokenType );

class Token : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~Token() override;
	virtual
	TokenType
	type() const = 0;

	virtual
	std::string
	to_string() const = 0;
};

class LiteralToken : public Token
{
public:
	LiteralToken();
	LiteralToken( numeric::Real value );

	
	TokenType
	type() const override;

	
	std::string
	to_string() const override;

	numeric::Real
	value() const;

	void
	value( numeric::Real setting );
private:
	numeric::Real value_;

};

class VariableToken : public Token
{
public:
	VariableToken();
	VariableToken( std::string const & name );

	
	TokenType
	type() const override;

	
	std::string
	to_string() const override;

	std::string name() const;
	void name( std::string const & name );
private:
	std::string name_;
};

class FunctionToken : public Token
{
public:
	FunctionToken();
	FunctionToken( std::string const & name, numeric::Size nargs );

	
	TokenType
	type() const override;

	
	std::string
	to_string() const override;

	std::string name() const;
	void name( std::string const & name );

	numeric::Size nargs() const;
	void nargs( Size setting);
private:
	std::string name_;
	numeric::Size nargs_;
};

class SimpleToken : public Token
{
public:
	SimpleToken();
	SimpleToken( TokenType type );

	
	TokenType
	type() const override;

	
	std::string
	to_string() const override;

	void
	set_token_type( TokenType type );

private:
	TokenType type_;
};

class TokenSet : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~TokenSet() override;
	TokenSet();

	void append( TokenCOP token );

	/// @brief Are there no more tokens remaining?
	bool empty() const;

	TokenCOP top() const;

	void pop();

	/// @brief Print contents after parser has encountered an error.
	void log_error() const;
	std::string print() const;
protected:

	/// @brief in the event of an error message, print the tokens up to the current token.
	void
	print_to_curr_pos() const;

private:
	std::list< TokenCOP > tokens_;
	std::list< TokenCOP >::iterator curr_pos_;
};

class ArithmeticScanner : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ArithmeticScanner() override;
	/// @brief Constructor which adds the "standard" set of min, max and sqrt functions.
	ArithmeticScanner();
	/// @brief Constructor which does not add the "standard" set of min, max and sqrt functions.
	ArithmeticScanner( bool );
	/// @brief Add the functions min, max and sqrt.
	void add_standard_functions();

	void add_variable( std::string const & name );
	void add_function( std::string const & name, numeric::Size nargs );
	TokenSetOP scan( std::string const & input_string );

private:
	LiteralTokenOP scan_literal( std::string const & input_string ) const;
	TokenOP scan_identifier( std::string const & input_string ) const;

	/// @brief print the contents of functions_ and variables_ to std error
	void log_error() const;

private:
	std::map< std::string, numeric::Size > functions_;
	std::map< std::string, numeric::Size > variables_; // Size is a placeholder; the name string is the only important data.

};

class BooleanExpressionScanner : public ArithmeticScanner
{
public:
	/// Simply adds more functions in its constructor so that
	/// a mixture of boolean and arithmetic expressions may be scanned.
	BooleanExpressionScanner();
};

///////////////////////////////////////////
//// BEGIN ABSTRAT SYNTAX TREE CLASSES ////
///////////////////////////////////////////


/// @brief Base class for Abstract Syntax Tree (AST) for the simple
/// Arithmetic language defined here.
class ArithmeticASTNode : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ArithmeticASTNode() override;
	virtual
	void
	visit( ASTVisitor & visitor ) const = 0;
};

class ArithmeticASTExpression : public ArithmeticASTNode
{
public:
	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_begin() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_end() const;

private:
	std::list< ArithmeticASTNodeCOP > children_;

};

class ArithmeticASTFunction : public ArithmeticASTNode
{
public:
	~ArithmeticASTFunction() override;

	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	FunctionTokenCOP
	function() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_begin() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_end() const;

private:
	FunctionTokenCOP function_;
	std::list< ArithmeticASTNodeCOP > children_;

};

class ArithmeticASTTerm : public ArithmeticASTNode
{
public:
	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_begin() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_end() const;

private:
	std::list< ArithmeticASTNodeCOP > children_;

};

class ArithmeticASTFactor : public ArithmeticASTNode
{
public:
	~ArithmeticASTFactor() override;

	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	ArithmeticASTNodeCOP
	child() const;
private:
	ArithmeticASTNodeCOP child_;

};

/// @brief either a variable or a literal.
class ArithmeticASTValue : public ArithmeticASTNode
{
public:
	ArithmeticASTValue();
	~ArithmeticASTValue() override;

	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	bool is_literal() const;
	numeric::Real literal_value() const;

	std::string variable_name() const;

private:
	bool is_literal_;
	numeric::Real literal_value_;
	std::string variable_name_;
};

class ArithmeticASTRestTerm : public ArithmeticASTNode
{
public:
	ArithmeticASTRestTerm();

	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_begin() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_end() const;

	TokenType
	rest_term_token() const;

private:
	TokenType rest_term_token_;
	std::list< ArithmeticASTNodeCOP > children_;

};

class ArithmeticASTRestExpression : public ArithmeticASTNode
{
public:
	ArithmeticASTRestExpression();

	
	void
	visit( ASTVisitor & visitor ) const override;

	virtual
	void
	parse( TokenSet & tokens );

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_begin() const;

	std::list< ArithmeticASTNodeCOP >::const_iterator
	children_end() const;

	TokenType
	rest_expression_token() const;

private:
	TokenType rest_expression_token_;
	std::list< ArithmeticASTNodeCOP > children_;

};

/// @brief Double-dispatch visitor pattern for abstract syntax tree
class ASTVisitor : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ASTVisitor() override;

	virtual
	void
	visit( ArithmeticASTExpression const & ) = 0;


	virtual
	void
	visit( ArithmeticASTFunction const & ) = 0;

	virtual
	void
	visit( ArithmeticASTTerm const & ) = 0;

	virtual
	void
	visit( ArithmeticASTFactor const & ) = 0;

	virtual
	void
	visit( ArithmeticASTValue const & ) = 0;

	virtual
	void
	visit( ArithmeticASTRestTerm const & ) = 0;

	virtual
	void
	visit( ArithmeticASTRestExpression const & ) = 0;

	virtual
	void
	visit( ArithmeticASTNode const & ) = 0;

};

/// Traverse the AST and print it to standard out
class ASTPrinter : public ASTVisitor
{
public:
	ASTPrinter();

	
	void
	visit( ArithmeticASTExpression const & ) override;


	
	void
	visit( ArithmeticASTFunction const & ) override;

	
	void
	visit( ArithmeticASTTerm const & ) override;

	
	void
	visit( ArithmeticASTFactor const & ) override;

	
	void
	visit( ArithmeticASTValue const & ) override;

	
	void
	visit( ArithmeticASTRestTerm const & ) override;

	
	void
	visit( ArithmeticASTRestExpression const & ) override;

	
	void
	visit( ArithmeticASTNode const & ) override;

	std::string
	ast_string( ArithmeticASTNode const & node );

	void
	pretty( bool setting );

private:
	void increment_indentation();
	void decrement_indentation();
	void indent_line();
	void finish_indented_line();

private:
	bool pretty_;
	Size indentation_level_;
	std::ostringstream ostrstream_;
	mutable bool last_dispatch_to_unknown_type_;

};

/// @brief Class to traverse the abstract syntax tree
/// produced by the parsing of a properly-formed string in the
/// Arithmetic expression language.  Produces an Expression tree
/// capable of performing arithmetic.  Connects the "variable" nodes
/// in this tree to the owning WrapperOptEMultifunc so that their
/// values can be retrieved during expression evaluation inside
/// the WrapperOptEMultifunc functor.
class ExpressionCreator : public ASTVisitor
{
public:
	ExpressionCreator();
	~ExpressionCreator() override;

	
	void
	visit( ArithmeticASTExpression const & ) override;

	
	void
	visit( ArithmeticASTFunction const & ) override;

	
	void
	visit( ArithmeticASTTerm const & ) override;

	
	void
	visit( ArithmeticASTFactor const & ) override;

	
	void
	visit( ArithmeticASTValue const & ) override;

	
	void
	visit( ArithmeticASTRestTerm const & ) override;

	
	void
	visit( ArithmeticASTRestExpression const & ) override;

	
	void
	visit( ArithmeticASTNode const & ) override;

	ExpressionCOP
	create_expression_tree( ArithmeticASTExpression const & );

	/// @brief Factory method to be implemented by derived classes
	/// which may wish to handle variable expressions in a specific manner
	virtual
	ExpressionCOP
	handle_variable_expression( ArithmeticASTValue const & );

	/// @brief Factory method to be implemented by derived classes
	/// which may wish to handle function expressions in a specific manner
	virtual
	ExpressionCOP
	handle_function_expression(
		FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	);

private:

	/// Created inside a traversal of the AST -- return expression trees
	/// using this variable, and keep living trees alive
	/// by storing them on the stack during the AST recursive traversal.
	ExpressionCOP last_constructed_expression_;
	ExpressionCOP semi_constructed_expression_;

};

class SimpleExpressionCreator : public ExpressionCreator
{
public:
	SimpleExpressionCreator();
	SimpleExpressionCreator( std::list< std::string > const & varnames );

	void
	add_variable( std::string const & varname );

	
	ExpressionCOP
	handle_variable_expression( ArithmeticASTValue const & ) override;

	VariableExpressionOP
	get_variable( std::string const & varname );

	std::map< std::string, VariableExpressionOP >
	variables() const;

private:
	std::map< std::string, VariableExpressionOP > variables_;

};

class BooleanExpressionCreator : public SimpleExpressionCreator
{
public:
	typedef SimpleExpressionCreator parent;
	typedef ExpressionCreator grandparent;
public:
	BooleanExpressionCreator();
	BooleanExpressionCreator( std::list< std::string > const & varnames );

	
	ExpressionCOP
	handle_function_expression(
		FunctionTokenCOP function,
		utility::vector1< ExpressionCOP > const & args
	) override;

};


///////////////////////////////////////////
///// END ABSTRAT SYNTAX TREE CLASSES /////
///////////////////////////////////////////

/// @brief Pure virtual base class to define arbitrary expressions for
/// scripting arithmetic operations (e.g. addition and multipliction).
class Expression : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~Expression() override;

	virtual
	numeric::Real
	operator() () const = 0;

	/// @brief Returns the expression for the partial derivative of this expression
	/// by the variable named varname.  If the partial derivative is always zero with respect
	/// to varname, returns null.
	virtual
	ExpressionCOP
	differentiate( std::string const & varname ) const = 0;

	virtual
	std::list< std::string >
	active_variables() const = 0;

};

class LiteralExpression : public Expression
{
public:
	LiteralExpression();
	LiteralExpression( numeric::Real value );

	void set_value( numeric::Real value );

	
	numeric::Real
	operator() () const override;

	/// @brief Returns null, since the derivative for a literal is always zero
	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

private:
	numeric::Real value_;

};

class VariableExpression : public Expression
{
public:
	VariableExpression( std::string const & name );
	VariableExpression( std::string const & name, numeric::Real value );

	
	numeric::Real
	operator() () const override;

	void set_value( numeric::Real value );

	std::string name() const;

	/// @brief Returns the literal expression 1 if name_ == varname_ and null otherwise
	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

private:
	std::string name_;
	numeric::Real value_;
};


class UnaryExpression : public Expression
{
public:
	UnaryExpression();
	~UnaryExpression() override;
	UnaryExpression( ExpressionCOP ex );

	void set_expression( ExpressionCOP ex );

	
	std::list< std::string >
	active_variables() const override;

public:
	ExpressionCOP ex() const;

private:
	ExpressionCOP ex_;

};

class BinaryExpression : public Expression
{
public:
	BinaryExpression();
	~BinaryExpression() override;
	BinaryExpression( ExpressionCOP e1, ExpressionCOP e2 );

	void set_first_expression( ExpressionCOP e1 );
	void set_second_expression( ExpressionCOP e2 );

	
	std::list< std::string >
	active_variables() const override;

public:
	ExpressionCOP e1() const;
	ExpressionCOP e2() const;

private:
	ExpressionCOP e1_;
	ExpressionCOP e2_;

};

class NaryExpression : public Expression
{
public:
	// undefinded NaryExpression();
	// undefinded NaryExpression( utility::vector1< ExpressionCOP > const & expressions );

	//void set_expressions( utility::vector1< ExpressionCOP > const & expressions );
	//Size size() const;
	//ExpressionCOP get_expression( Size id ) const;

private:
	//Size nargs_;
	utility::vector1< ExpressionCOP > expressions_;
};

class SquarerootExpression : public UnaryExpression
{
public:
	SquarerootExpression();
	SquarerootExpression( ExpressionCOP ex );

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class AbsoluteValueExpression : public UnaryExpression
{
public:
	AbsoluteValueExpression();
	AbsoluteValueExpression( ExpressionCOP ex );

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class AddExpression : public BinaryExpression
{
public:
	AddExpression();
	AddExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the sum of expression 1 and expression 2.
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class SubtractExpression : public BinaryExpression
{
public:
	SubtractExpression();
	SubtractExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the difference between expression 1 and expression 2.
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};


class MultiplyExpression : public BinaryExpression
{
public:
	MultiplyExpression();
	MultiplyExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the product of expression 1 and expression 2
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class DivideExpression : public BinaryExpression
{
public:
	DivideExpression();
	DivideExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the quotient of expression 1 and expression 2
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class MaxExpression : public BinaryExpression
{
public:
	MaxExpression();
	MaxExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the max of e1 and e2.
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

};

class MinExpression : public BinaryExpression
{
public:
	MinExpression();
	MinExpression( ExpressionCOP e1, ExpressionCOP e2 );

	/// @brief Returns the min of e1 and e2
	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

};


/// @brief Evaluates ee1 when e1 is larger than e2; evaluates ee2 otherwise.
class MetaMaxExpression : public Expression
{
public:
	//MetaMaxExpression();
	MetaMaxExpression( ExpressionCOP e1, ExpressionCOP e2, ExpressionCOP ee1, ExpressionCOP ee2 );

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

private:
	ExpressionCOP e1_, e2_;
	ExpressionCOP ee1_, ee2_;

};

/// @brief Evaluates ee1 when e1 is less than e2; evaluates ee2 otherwise.
class MetaMinExpression : public Expression
{
public:
	//MetaMinExpression();
	MetaMinExpression( ExpressionCOP e1, ExpressionCOP e2, ExpressionCOP ee1, ExpressionCOP ee2 );

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	
	std::list< std::string >
	active_variables() const override;

private:
	ExpressionCOP e1_, e2_;
	ExpressionCOP ee1_, ee2_;

};

/// BEGIN PSEUDO-BOOLEAN EXPRESSIONS
/// DO NOT TRY TO DIFFERENTIATE THEM
/// 1. Arithmetic Comparison Operators
class EqualsExpression : public BinaryExpression
{
public:
	EqualsExpression();
	EqualsExpression( ExpressionCOP e1, ExpressionCOP e2 );
	~EqualsExpression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

/// Greater Than
class GT_Expression : public BinaryExpression
{
public:
	GT_Expression();
	GT_Expression( ExpressionCOP e1, ExpressionCOP e2 );
	~GT_Expression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

/// Greater Than or Equal To
class GTE_Expression : public BinaryExpression
{
public:
	GTE_Expression();
	GTE_Expression( ExpressionCOP e1, ExpressionCOP e2 );
	~GTE_Expression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

/// Less Than
class LT_Expression : public BinaryExpression
{
public:
	LT_Expression();
	LT_Expression( ExpressionCOP e1, ExpressionCOP e2 );
	~LT_Expression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

/// Less Than or Equal To
class LTE_Expression : public BinaryExpression
{
public:
	LTE_Expression();
	LTE_Expression( ExpressionCOP e1, ExpressionCOP e2 );
	~LTE_Expression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};


/// 2. Boolean Logic Operators
class AndExpression : public BinaryExpression
{
public:
	AndExpression();
	AndExpression( ExpressionCOP e1, ExpressionCOP e2 );
	~AndExpression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class OrExpression : public BinaryExpression
{
public:
	OrExpression();
	OrExpression( ExpressionCOP e1, ExpressionCOP e2 );
	~OrExpression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class NotExpression : public UnaryExpression
{
public:
	NotExpression();
	NotExpression( ExpressionCOP e );
	~NotExpression() override;

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

};

class ITEExpression : public Expression
{
public:

	ITEExpression(
		ExpressionCOP condition,
		ExpressionCOP then_clause,
		ExpressionCOP else_clause
	);

	
	numeric::Real
	operator() () const override;

	
	ExpressionCOP
	differentiate( std::string const & varname ) const override;

	ExpressionCOP
	condition() const;

	ExpressionCOP
	then_expression() const;

	ExpressionCOP
	else_expression() const;

	
	std::list< std::string >
	active_variables() const override;


private:
	ExpressionCOP condition_;
	ExpressionCOP then_expression_;
	ExpressionCOP else_expression_;

};


ExpressionCOP
parse_string_to_expression( std::string const & input_string );

ExpressionCOP
parse_string_to_boolean_expression( std::string const & input_string );


}
}

#endif
