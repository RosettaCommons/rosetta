// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/expression_parser/Arithmetic.cc
/// @brief  Parse tree for arithmetic operations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <numeric/expression_parser/Arithmetic.hh>

/// Core headers
#include <numeric/types.hh>

/// Utility headers
#include <utility>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

/// C++ headers
#include <iostream>
#include <cmath>

#include <utility/vector1.hh>


namespace numeric {
namespace expression_parser {

/// @details Auto-generated virtual destructor
Expression::~Expression() = default;

/// @details Auto-generated virtual destructor
ASTVisitor::~ASTVisitor() = default;

/// @details Auto-generated virtual destructor
ArithmeticASTNode::~ArithmeticASTNode() = default;

/// @details Auto-generated virtual destructor
ArithmeticScanner::~ArithmeticScanner() = default;

/// @details Auto-generated virtual destructor
TokenSet::~TokenSet() = default;

/// @details Auto-generated virtual destructor
Token::~Token() = default;

std::string
token_type_name( TokenType tt )
{
	switch ( tt ) {
	case INVALID_TOKEN_TYPE :
		return "INVALID_TOKEN_TYPE";
	case LITERAL :
		return "LITERAL";
	case VARIABLE :
		return "VARIABLE";
	case FUNCTION :
		return "FUNCTION";
	case COMMA :
		return "COMMA";
	case LEFT_PAREN :
		return "LEFT_PAREN";
	case RIGHT_PAREN :
		return "RIGHT_PAREN";
	case PLUS_SYMBOL :
		return "PLUS_SYMBOL";
	case SUBTRACT_SYMBOL :
		return "SUBTRACT_SYMBOL";
	case MULTIPLY_SYMBOL :
		return "MULTIPLY_SYMBOL";
	case DIVIDE_SYMBOL :
		return "DIVIDE_SYMBOL";
	default :
		return "ERROR IN token_type_name -- unrecognized token type!";
	}

	return "";
}

bool
is_numeral( char c ) {
	return c >= '0' && c <= '9';
}

bool
is_letter( char c ) {
	return ( c >= 'a' && c <= 'z' ) || ( c >= 'A' && c <= 'Z' );
}

LiteralToken::LiteralToken() : value_( 0.0 ) {}


LiteralToken::LiteralToken( numeric::Real value ) : value_( value ) {}

TokenType
LiteralToken::type() const
{
	return LITERAL;
}

std::string
LiteralToken::to_string() const
{
	return "LITERAL(" + utility::to_string( value_ ) + ")";
}

numeric::Real
LiteralToken::value() const
{
	return value_;
}

void
LiteralToken::value( numeric::Real setting )
{
	value_ = setting;
}


VariableToken::VariableToken() : name_() {}
VariableToken::VariableToken( std::string  name ) : name_(std::move( name )) {}

TokenType
VariableToken::type() const
{
	return VARIABLE;
}

std::string
VariableToken::to_string() const
{
	return "VARIABLE(" + name_ + ")";
}

std::string
VariableToken::name() const
{
	return name_;
}

void VariableToken::name( std::string const & name )
{
	name_ = name;
}

FunctionToken::FunctionToken() : name_(), nargs_( 0 ) {}
FunctionToken::FunctionToken( std::string  name, numeric::Size nargs ) :
	name_(std::move( name )),
	nargs_( nargs )
{}

TokenType
FunctionToken::type() const
{
	return FUNCTION;
}

std::string
FunctionToken::to_string() const
{
	return "FUNCTION(" + name_ + "," + utility::to_string(nargs_) + ")";
}

std::string
FunctionToken::name() const
{
	return name_;
}

void
FunctionToken::name( std::string const & setting )
{
	name_ = setting;
}

numeric::Size
FunctionToken::nargs() const
{
	return nargs_;
}

void FunctionToken::nargs( Size setting ) {
	nargs_ = setting;
}


SimpleToken::SimpleToken() : type_( INVALID_TOKEN_TYPE ) {}

SimpleToken::SimpleToken( TokenType type ) : type_( type ) {}

TokenType
SimpleToken::type() const
{
	return type_;
}

std::string
SimpleToken::to_string() const
{
	return token_type_name( type_ );
}


void
SimpleToken::set_token_type( TokenType setting )
{
	type_ = setting;
}


TokenSet::TokenSet() = default;

/// @details Appends the new token to the list and resets
/// the token iterator to the beginning of the list
void TokenSet::append( TokenCOP token )
{
	tokens_.push_back( token );
	curr_pos_ = tokens_.begin();
}

bool
TokenSet::empty() const {
	return tokens_.empty() || curr_pos_ == tokens_.end();
}

TokenCOP
TokenSet::top() const
{
	if ( curr_pos_ == tokens_.end() ) {
		utility_exit_with_message("No tokens remaining" );
	}
	return *curr_pos_;
}

void TokenSet::pop()
{
	if ( curr_pos_ == tokens_.end() ) {
		print_to_curr_pos();
		utility_exit_with_message("Invalid token pop" );
	}
	++curr_pos_;
}

void
TokenSet::log_error() const
{
	print_to_curr_pos();
}

std::string
TokenSet::print() const
{
	std::ostringstream ostream;
	ostream << "Tokens:\n";
	for ( auto const & token : tokens_ ) {
		ostream << token->to_string() << "\n";
	}
	return ostream.str();
}

void
TokenSet::print_to_curr_pos() const
{
	std::cerr << "TokenSet parsing error" << std::endl;
	Size count( 0 );
	for ( std::list< TokenCOP >::const_iterator
			iter = tokens_.begin(), iter_end = curr_pos_;
			iter != iter_end; ++iter ) {
		++count;
		std::cerr << count << ": " << (*iter)->to_string() << std::endl;
	}

}

ArithmeticScanner::ArithmeticScanner()
{
	add_standard_functions();
}

ArithmeticScanner::ArithmeticScanner( bool /*dummy */ )
{
	/// no op -- do not add the standard functions.
}

void ArithmeticScanner::add_standard_functions()
{
	/// Three functions "built in."  min, max and sqrt.
	add_function( "min", 2 );
	add_function( "max", 2 );
	add_function( "sqrt", 1 );
}

void ArithmeticScanner::add_variable( std::string const & name )
{
	if ( name.size() == 0 ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_variable.  Illegal variable name with length 0" );
	} else if ( is_numeral( name[ 0 ] ) ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_variable.  Illegal variable name beginning with a numeral: " + name );
	}
	if ( variables_.find( name ) != variables_.end() ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_variable.  Variable name has already been added: " + name );
	}
	if ( functions_.find( name ) != functions_.end() ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_variable.  Variable name conflicts with already added function name: " + name );
	}
	variables_.insert( std::make_pair( name, 1 ));
}

void ArithmeticScanner::add_function(
	std::string const & name,
	numeric::Size nargs
)
{
	if ( name.size() == 0 ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_function.  Illegal function name with length 0" );
	} else if ( is_numeral( name[ 0 ] ) ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_function.  Illegal function name beginning with a numeral: " + name );
	} else if ( nargs == 0 ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_function.  Illegal function with no arguments: " + name );
	}
	if ( variables_.find( name ) != variables_.end() ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_function.  Function name has already conflicts with already added variable name: " + name );
	}
	if ( functions_.find( name ) != functions_.end() ) {
		utility_exit_with_message( "Error in ArithmeticScanner::add_variable.  Function name has already been added: " + name );
	}

	functions_.insert( std::make_pair( name,  nargs ));
}

TokenSetOP
ArithmeticScanner::scan( std::string const & input_string )
{
	TokenSetOP tokens( new TokenSet );
	Size pos_token_begin( 0 ), pos_curr( 0 );
	bool scanning_literal = false;
	while ( pos_curr <= input_string.size() ) {

		if ( pos_token_begin != pos_curr ) {
			/// we're in the middle of a token
			if ( scanning_literal ) {
				bool done = false;
				if ( pos_curr == input_string.size() ) {
					done = true;
				} else {
					char curr_char = input_string[ pos_curr ];
					if ( curr_char == ' ' ||  curr_char == '(' ||
							curr_char == ')' || curr_char == ',' ||
							curr_char == '+' || curr_char == '*' ||
							curr_char == '/' || curr_char == '\t' ) {
						done = true;
					} else if ( curr_char == '-' ) {
						done = true;
						/// exceptions! e.g. 1e-6
						assert( pos_curr > 0 ); // so that subtracting 1 doesn't wrap to Size(-1) (4 billion on 32-bit machines).
						if ( pos_token_begin + 1 != pos_curr &&
								input_string[ pos_curr - 1 ] == 'e' &&
								pos_curr + 1 < input_string.size() &&
								is_numeral( input_string[ pos_curr + 1 ] ) ) {
							done = false;
						}
					}
				}

				if ( done ) {
					std::string literal_string = input_string.substr( pos_token_begin, pos_curr - pos_token_begin );
					tokens->append( scan_literal( literal_string ) );
					pos_token_begin = pos_curr;
					--pos_curr; // incremented back to its current value
				}

			} else { // scanning a variable or a function
				bool done = false;
				if ( pos_curr == input_string.size() ) {
					done = true;
				} else {
					char curr_char = input_string[ pos_curr ];
					if ( curr_char == ' ' ||  curr_char == '(' ||
							curr_char == ')' || curr_char == ',' ||
							curr_char == '+' || curr_char == '*' ||
							curr_char == '/' || curr_char == '-' || curr_char == '\t' ) {
						done = true;
					}
				}

				if ( done ) {
					std::string identifier_string = input_string.substr( pos_token_begin, pos_curr - pos_token_begin );
					tokens->append( scan_identifier( identifier_string ) );
					pos_token_begin = pos_curr;
					--pos_curr; // incremented back to its current value
				}
			}

		} else if ( pos_curr == input_string.size() ) {
			break; /// Finished tokenizing the input string and we're not in the middle of an untokenized literal or identifier; quit.
		} else if ( input_string[ pos_curr ] != ' ' && input_string[ pos_curr ] != '\t' ) {
			scanning_literal = false;

			if ( input_string[ pos_curr ] == '(' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( LEFT_PAREN ) ) ) );
				++pos_token_begin;
			} else if ( input_string[ pos_curr ] == ')' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( RIGHT_PAREN ) ) ) );
				++pos_token_begin;
			} else if ( input_string[ pos_curr ] == ',' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( COMMA ) ) ) );
				++pos_token_begin;
			} else if ( input_string[ pos_curr ] == '+' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( PLUS_SYMBOL ) ) ) );
				++pos_token_begin;
			} else if ( input_string[ pos_curr ] == '-' ) {
				if ( pos_curr + 1 < input_string.size() &&
						(is_numeral(input_string[ pos_curr + 1 ]) || input_string[ pos_curr + 1 ] == '.') ) {
					scanning_literal = true;
					pos_token_begin = pos_curr;
				} else {
					tokens->append( TokenCOP( TokenOP( new SimpleToken( SUBTRACT_SYMBOL ) ) ));
					++pos_token_begin;
				}
			} else if ( input_string[ pos_curr ] == '*' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( MULTIPLY_SYMBOL ) ) ) );
				++pos_token_begin;
			} else if ( input_string[ pos_curr ] == '/' ) {
				tokens->append( TokenCOP( TokenOP( new SimpleToken( DIVIDE_SYMBOL ) ) ) );
				++pos_token_begin;
			} else if ( is_numeral( input_string[ pos_curr ] ) ) {
				scanning_literal = true;
				pos_token_begin = pos_curr;
			} else if ( is_letter( input_string[ pos_curr ] ) || input_string[ pos_curr ] == '$' ) {
				pos_token_begin = pos_curr;
			} else {
				utility_exit_with_message( "Illegal character encountered in scanning the following string near position "
					+ utility::to_string( pos_curr ) + " = '" + utility::to_string( input_string[ pos_curr ] ) +  "'\n" + input_string );
			}
		} else {
			/// input character is a ' ' or a '\t' and we aren't yet scanning a literal or identifier
			/// increment pos_token_begin
			++pos_token_begin;
		}
		++pos_curr;

	}

	return tokens;
}


LiteralTokenOP
ArithmeticScanner::scan_literal( std::string const & input_string ) const
{
	bool found_e = false; bool found_point = false;
	for ( Size ii = 0; ii < input_string.size(); ++ii ) {
		if ( ! is_numeral( input_string[ ii ] ) ) {
			if ( input_string[ ii ] == 'e' || input_string[ ii ] == 'E' ) {
				if ( ! found_e ) {
					found_e = true;
					continue;
				} else {
					utility_exit_with_message( "Error in trying to parse the following string as a numeric literal -- too many letter e's\n"
						+ utility::to_string( ii ) + " " + utility::to_string( input_string[ ii ] ) + "\n" + input_string );
				}
			}
			if ( input_string[ ii ] == '.' ) {
				if ( ! found_point ) {
					found_point = true;
					continue;
				} else {
					utility_exit_with_message( "Error in trying to parse the following string as a numeric literal -- too many periods\n"
						+ utility::to_string( ii ) + " " + utility::to_string( input_string[ ii ] ) + "\n" + input_string );
				}
			}
			if ( input_string[ ii ] == '-' ) {
				if ( ii == 0 ||
						( ii + 1 < input_string.size() && is_numeral( input_string[ ii + 1 ] ) &&
						( input_string[ ii-1 ] == 'e' || input_string[ ii-1 ] == 'E' )) ) {
					continue;
				}
			}
			utility_exit_with_message( "Error trying to parse the following string as a numeric literal at position "
				+ utility::to_string( ii ) + " = '" + utility::to_string( input_string[ ii ] ) + "'\n" + input_string );
		}
	}
	numeric::Real value = utility::from_string( input_string, numeric::Real( 0.0 ) );
	return LiteralTokenOP( new LiteralToken( value ) );
}

TokenOP
ArithmeticScanner::scan_identifier( std::string const & input_string ) const
{
	if ( input_string.size() == 0 ) {
		utility_exit_with_message( "Error in scan_identifier: input_string is empty" );
	}

	if ( ! is_letter( input_string[ 0 ] ) && input_string[ 0 ] != '$' ) {
		utility_exit_with_message( "Error in scan_identifier: input_string must begin with a letter or a '$' character: '" + utility::to_string( input_string[ 0 ] ) +"'/n" + input_string  );
	}

	for ( Size ii = 1; ii < input_string.size(); ++ii ) {
		char iichar = input_string[ ii ];
		if ( ! is_numeral( iichar ) && ! is_letter( iichar ) && iichar != '_' && iichar != '$' ) {
			utility_exit_with_message( "Error trying to parse the following string as an identifier at position "
				+ utility::to_string( ii ) + " = '" + utility::to_string( input_string[ ii ] ) + "'\n" + input_string );
		}
	}

	/// try inserting as a variable
	if ( variables_.find( input_string ) != variables_.end() ) {
		return TokenOP( new VariableToken( input_string ) );
	}
	/// ok -- now try inserting as a function
	if ( functions_.find( input_string ) != functions_.end() ) {
		return TokenOP( new FunctionToken( input_string, functions_.find( input_string )->second ) );
	}

	log_error();
	utility_exit_with_message( "Error in scan_identifier: Failed to find input_string in either variables_ or functions_ maps\n" + input_string );
	return nullptr; /// appease compiler
}

/// @brief print the contents of functions_ and variables_ to std error
void ArithmeticScanner::log_error() const
{
	std::cerr << "An error has occurred.  ArithmeticScanner contents: " << std::endl;
	for ( auto const & variable : variables_ ) {
		std::cerr << "variable: " << variable.first << std::endl;
	}
	for ( auto const & function : functions_ ) {
		std::cerr << "functions: " << function.first << " #args: " << function.second << std::endl;
	}
	std::cerr << std::endl;
}


BooleanExpressionScanner::BooleanExpressionScanner() : ArithmeticScanner()
{
	add_function( "EQUALS", 2 );
	add_function( "GT", 2 );
	add_function( "GTE", 2 );
	add_function( "LT", 2 );
	add_function( "LTE", 2 );
	add_function( "AND", 2 );
	add_function( "OR", 2 );
	add_function( "NOT", 1 );
}

///////////////////////////////////////////
//// BEGIN ABSTRAT SYNTAX TREE CLASSES ////
///////////////////////////////////////////

/// EXPRESSION

void
ArithmeticASTExpression::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTExpression::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) {
		tokens.log_error();
		utility_exit_with_message( "Error parsing ArithmeticASTExpression.  Empty token list!" );
	}

	//std::cout << "ASTExpression " << token_type_name( tokens.top()->type() );
	ArithmeticASTTermOP child1( new ArithmeticASTTerm );
	child1->parse( tokens );
	children_.push_back( child1 );
	ArithmeticASTRestExpressionOP child2( new ArithmeticASTRestExpression );
	child2->parse( tokens );
	children_.push_back( child2 );
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTExpression::children_begin() const
{
	return children_.begin();
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTExpression::children_end() const
{
	return children_.end();
}

//// FUNCTION

ArithmeticASTFunction::~ArithmeticASTFunction() = default;

void
ArithmeticASTFunction::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTFunction::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) {
		tokens.log_error();
		utility_exit_with_message( "Error parsing ArithmeticASTFunction.  Empty token list!" );
	}

	if ( ! tokens.empty() ) {
		//std::cout << "ArithmeticASTFunction " << token_type_name( tokens.top()->type() );
		if ( tokens.top()->type() != FUNCTION ) {
			tokens.log_error();
			utility_exit_with_message( "Error in ArithmeticASTFunction parse: non function input!" );
		}
		TokenCOP toptoken = tokens.top();
		function_ = FunctionTokenCOP( utility::pointer::dynamic_pointer_cast< numeric::expression_parser::FunctionToken const > ( toptoken ) );
		if ( ! function_ ) {
			utility_exit_with_message( "Error.  Dynamic cast from token claiming to be FUNCTION type failed" );
		}
		Size func_nargs = function_->nargs();
		if ( func_nargs == 0 ) {
			utility_exit_with_message("Error in ArithmeticASTFunction parse: function takes no arguments" );
		}
		assert( tokens.top()->type() == FUNCTION );
		tokens.pop();

		if ( tokens.top()->type() != LEFT_PAREN ) {
			tokens.log_error();
			utility_exit_with_message( "Error parsing ArithmeticASTFunction: expected left paren.");
		}
		assert( tokens.top()->type() == LEFT_PAREN );
		tokens.pop();

		for ( Size ii = 1; ii < func_nargs; ++ii ) {
			ArithmeticASTExpressionOP exp( new ArithmeticASTExpression );
			exp->parse( tokens );
			if ( tokens.top()->type() != COMMA ) {
				tokens.log_error();
				utility_exit_with_message( "Error in parsing ArithmeticASTFunction: expected comma." );
			}

			assert( tokens.top()->type() == COMMA );
			tokens.pop();

			children_.push_back( exp );
		}
		ArithmeticASTExpressionOP exp( new ArithmeticASTExpression );
		exp->parse( tokens );
		children_.push_back( exp );
		if ( tokens.top()->type() != RIGHT_PAREN ) {
			tokens.log_error();
			utility_exit_with_message( "Error in parsing ArithmeticASTFunction: expected right paren.");
		}
		tokens.pop();
	}
}

FunctionTokenCOP
ArithmeticASTFunction::function() const
{
	if ( function_ == nullptr ) {
		utility_exit_with_message( "function() request of ArithmeticASTFunction would return null pointer" );
	}
	return function_;
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTFunction::children_begin() const
{
	return children_.begin();
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTFunction::children_end() const
{
	return children_.end();
}

void
ArithmeticASTTerm::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}


void
ArithmeticASTTerm::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) {
		tokens.log_error();
		utility_exit_with_message("Error parsing ArithmeticASTTerm.  Empty token list!" );
	}
	ArithmeticASTFactorOP factor( new ArithmeticASTFactor );
	factor->parse( tokens );
	children_.push_back( factor );
	ArithmeticASTRestTermOP rest_term( new ArithmeticASTRestTerm );
	rest_term->parse( tokens );
	children_.push_back( rest_term );
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTTerm::children_begin() const
{
	return children_.begin();
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTTerm::children_end() const
{
	return children_.end();
}

ArithmeticASTFactor::~ArithmeticASTFactor() = default;

void
ArithmeticASTFactor::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTFactor::parse( TokenSet & tokens )
{
	TokenType toptype = tokens.top()->type();
	if ( toptype == FUNCTION ) {
		ArithmeticASTFunctionOP func( new ArithmeticASTFunction );
		func->parse( tokens );
		child_ = func;
	} else if ( toptype == VARIABLE || toptype == LITERAL ) {
		ArithmeticASTValueOP val( new ArithmeticASTValue );
		val->parse( tokens );
		child_ = val;
	} else if ( toptype == LEFT_PAREN ) {
		ArithmeticASTExpressionOP exp( new ArithmeticASTExpression );
		tokens.pop();
		exp->parse( tokens );
		child_ = exp;
		if ( tokens.top()->type() != RIGHT_PAREN ) {
			tokens.log_error();
			utility_exit_with_message( "Parse error in ArithmeticASTFactor, expected RIGHT_PAREN but got " + token_type_name( tokens.top()->type() ) );
		}
		tokens.pop();
	} else {
		tokens.log_error();
		utility_exit_with_message( "Parse error in ArithmeticASTFactor, expected a function, variable, literal or left-paren but got " + token_type_name( toptype ) );
	}
}

ArithmeticASTNodeCOP
ArithmeticASTFactor::child() const
{
	return child_;
}

ArithmeticASTValue::ArithmeticASTValue()
:
	is_literal_( false ),
	literal_value_( 0.0 ),
	variable_name_( "UNASSIGNED_VARIABLE_NAME" )
{}

ArithmeticASTValue::~ArithmeticASTValue() = default;

void
ArithmeticASTValue::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTValue::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) {
		tokens.log_error();
		utility_exit_with_message( "Error parsing ArithmeticASTValue.  Token set empty" );
	}
	TokenType toptoken = tokens.top()->type();
	if ( toptoken == VARIABLE ) {
		is_literal_ = false;
		VariableTokenCOP vartok = VariableTokenCOP( utility::pointer::dynamic_pointer_cast< numeric::expression_parser::VariableToken const > ( tokens.top() ) );
		if ( vartok->name() == "UNASSIGNED_VARIABLE_NAME" ) {
			utility_exit_with_message( "Illegal variable name in ArithmeticASTValue -- UNASSIGNED_VARIABLE_NAME is a reserved string" );
		}
		variable_name_ = vartok->name();
		tokens.pop();
	} else if ( toptoken == LITERAL ) {
		LiteralTokenCOP littok = LiteralTokenCOP( utility::pointer::dynamic_pointer_cast< numeric::expression_parser::LiteralToken const > ( tokens.top() ) );
		literal_value_ = littok->value();
		is_literal_ = true;
		tokens.pop();
	} else {
		tokens.log_error();
		utility_exit_with_message( "Error in parsing ArithmeticASTValue.  Expected variable or literal but got " + token_type_name( toptoken ) );
	}

}

bool ArithmeticASTValue::is_literal() const
{
	return is_literal_;
}

numeric::Real ArithmeticASTValue::literal_value() const
{
	if ( ! is_literal_ ) {
		utility_exit_with_message( "Illegal access of literal_value() for ArithmeticASTValue that is not a literal (and is instead a variable).");
	}
	return literal_value_;
}

std::string ArithmeticASTValue::variable_name() const
{
	if ( is_literal_ ) {
		utility_exit_with_message( "Illegal access of variable_name for ArithmeticASTValue that is not a variable (and is instead a literal).");
	}
	if ( variable_name_ == "UNASSIGNED_VARIABLE_NAME" ) {
		utility_exit_with_message( "Illegal access of variable_name for ArithmeticASTValue; variable name has not been assigned");
	}
	return variable_name_;
}

ArithmeticASTRestTerm::ArithmeticASTRestTerm() : rest_term_token_( INVALID_TOKEN_TYPE ) {}

void
ArithmeticASTRestTerm::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTRestTerm::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) {
		return;
	}
	TokenType toptoken = tokens.top()->type();
	if ( toptoken == MULTIPLY_SYMBOL || toptoken == DIVIDE_SYMBOL ) {
		rest_term_token_ = toptoken;
		tokens.pop();
		ArithmeticASTFactorOP factor( new ArithmeticASTFactor );
		factor->parse( tokens );
		children_.push_back( factor );
		ArithmeticASTRestTermOP restterm( new ArithmeticASTRestTerm );
		restterm->parse( tokens );
		children_.push_back( restterm );
	}
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTRestTerm::children_begin() const
{
	return children_.begin();
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTRestTerm::children_end() const
{
	return children_.end();
}

TokenType
ArithmeticASTRestTerm::rest_term_token() const
{
	if ( rest_term_token_ == INVALID_TOKEN_TYPE ) {
		/// Something is wrong -- is this an empty RestTerm?
		if ( children_.begin() != children_.end() ) {
			utility_exit_with_message( "Internal error in ArithmeticASTRestTerm: Invalid multiply/divide symbol, yet non-empty child list");
		} else {
			/// If so, then the calling function should not have tried to ask for the token type.
			utility_exit_with_message( "Error accessinig ArithmeticASTRestTerm rest_term_token; should not be accessed for childless object.");
		}
	}
	return rest_term_token_;
}

ArithmeticASTRestExpression::ArithmeticASTRestExpression()
:
	rest_expression_token_( INVALID_TOKEN_TYPE )
{}

void
ArithmeticASTRestExpression::visit( ASTVisitor & visitor ) const
{
	visitor.visit( *this );
}

void
ArithmeticASTRestExpression::parse( TokenSet & tokens )
{
	if ( tokens.empty() ) return;
	TokenType toptoken = tokens.top()->type();
	if ( toptoken == PLUS_SYMBOL || toptoken == SUBTRACT_SYMBOL ) {
		rest_expression_token_ = toptoken;
		tokens.pop();
		ArithmeticASTTermOP term( new ArithmeticASTTerm );
		term->parse( tokens );
		children_.push_back( term );
		ArithmeticASTRestExpressionOP restexpr( new ArithmeticASTRestExpression );
		restexpr->parse( tokens );
		children_.push_back( restexpr );
	}
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTRestExpression::children_begin() const
{
	return children_.begin();
}

std::list< ArithmeticASTNodeCOP >::const_iterator
ArithmeticASTRestExpression::children_end() const
{
	return children_.end();
}

TokenType
ArithmeticASTRestExpression::rest_expression_token() const
{
	if ( rest_expression_token_ == INVALID_TOKEN_TYPE ) {
		/// Something is wrong
		if ( children_.begin() != children_.end() ) {
			utility_exit_with_message( "Internal error in ArithmeticASTRestExpression: Invalid multiply/divide symbol, yet non-empty child list");
		} else {
			utility_exit_with_message( "Error accessinig ArithmeticASTRestExpression rest_expression_token; should not be accessed for childless object.");
		}
	}
	return rest_expression_token_;
}

///// END AST /////


///// BEGIN VISITORS /////
/// Traverse the AST and print it to standard out

ASTPrinter::ASTPrinter() :
	pretty_( true ),
	indentation_level_( 0 ),
	last_dispatch_to_unknown_type_( false )
{}

void
ASTPrinter::visit( ArithmeticASTExpression const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTExpression(";
	finish_indented_line();
	increment_indentation();
	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		finish_indented_line();
	}
	decrement_indentation();
	indent_line();
	ostrstream_ << ")";
}

void
ASTPrinter::visit( ArithmeticASTFunction const & node )
{
	last_dispatch_to_unknown_type_ = false;
	Size count_args( 0 );
	indent_line();
	ostrstream_ << "ArithmeticASTFunction:" << node.function()->name() << ":" << node.function()->nargs() << "(";
	finish_indented_line();
	increment_indentation();
	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		++count_args;
		(*iter)->visit( *this );
		if ( count_args != node.function()->nargs() ) {
			ostrstream_ << ",";
		}
		finish_indented_line();
	}
	decrement_indentation();
	indent_line();
	ostrstream_ << ")";
}

void
ASTPrinter::visit( ArithmeticASTTerm const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTTerm(";
	finish_indented_line();
	increment_indentation();
	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		finish_indented_line();
	}
	decrement_indentation();
	indent_line();
	ostrstream_ << ")";
}

void
ASTPrinter::visit( ArithmeticASTFactor const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTFactor( ";
	finish_indented_line();
	increment_indentation();
	node.child()->visit( *this );
	finish_indented_line();
	decrement_indentation();
	indent_line();
	ostrstream_ << ") ";
}

void
ASTPrinter::visit( ArithmeticASTValue const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTValue:";
	if ( node.is_literal() ) {
		ostrstream_ << "literal:" << node.literal_value();
	} else {
		ostrstream_ << "variable:" << node.variable_name();
	}
}

void
ASTPrinter::visit( ArithmeticASTRestTerm const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTRestTerm";
	if ( node.children_begin() == node.children_end() ) {
		ostrstream_ << "()";
	} else {
		ostrstream_ << ":" << token_type_name( node.rest_term_token() ) << "(";
		finish_indented_line();
		increment_indentation();
		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			(*iter)->visit( *this );
			finish_indented_line();
		}
		decrement_indentation();
		indent_line();
		ostrstream_ << ")";
	}
}

void
ASTPrinter::visit( ArithmeticASTRestExpression const & node )
{
	last_dispatch_to_unknown_type_ = false;
	indent_line();
	ostrstream_ << "ArithmeticASTRestExpression";
	if ( node.children_begin() == node.children_end() ) {
		ostrstream_ << "()";
	} else {
		ostrstream_ << ":" << token_type_name( node.rest_expression_token() ) << "(";
		finish_indented_line();
		increment_indentation();
		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			(*iter)->visit( *this );
			finish_indented_line();
		}
		decrement_indentation();
		indent_line();
		ostrstream_ << ")";
	}
}

void
ASTPrinter::visit( ArithmeticASTNode const & node )
{
	/// Danger, if we have an unkown type and dispatch to it and it dispatches back to this
	/// we enter an infinite recursion.  Quit as soon as dispatch fails twice in a row.
	if (  last_dispatch_to_unknown_type_ ) {
		std::cerr << "ERROR -- Infinite loop catch.  Two consequitive dispatch failures from ASTPrinter.  ASTNode derived class unhandled in dispatch!" << std::endl;
		utility_exit_with_message("Stuck in (likely) infinite loop");
	}
	last_dispatch_to_unknown_type_ = false;

	node.visit( *this );
}

std::string
ASTPrinter::ast_string( ArithmeticASTNode const & node )
{
	ostrstream_.str( "" );
	node.visit( *this );
	return ostrstream_.str();
}

void
ASTPrinter::pretty( bool setting )
{
	pretty_ = setting;
}


void ASTPrinter::increment_indentation()
{
	++indentation_level_;
}

void ASTPrinter::decrement_indentation()
{
	--indentation_level_;
}

void ASTPrinter::indent_line()
{
	if ( pretty_ ) {
		for ( Size ii = 1; ii <= indentation_level_; ++ii ) {
			ostrstream_ << "  ";
		}
	}
}

void ASTPrinter::finish_indented_line()
{
	if ( pretty_ ) {
		ostrstream_ << "\n";
	} else {
		ostrstream_ << " ";
	}
}

//// Expression visitor
ExpressionCreator::ExpressionCreator() = default;

ExpressionCreator::~ExpressionCreator() = default;

void
ExpressionCreator::visit( ArithmeticASTExpression const & node)
{
	utility::vector1< ExpressionCOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		if ( last_constructed_expression_ ) {
			expressions.push_back( last_constructed_expression_ );
			semi_constructed_expression_ = last_constructed_expression_;
		}
	}

	if ( expressions.size() == 1 ) {
		last_constructed_expression_ = expressions[ 1 ];
		semi_constructed_expression_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_expression_ = expressions[ 2 ];
		semi_constructed_expression_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTExpression.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}

}

void
ExpressionCreator::visit( ArithmeticASTFunction const & node)
{
	utility::vector1< ExpressionCOP > args( node.function()->nargs(), nullptr );
	Size count_args = 0;
	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		++count_args;
		(*iter)->visit( *this );
		args[ count_args ] = last_constructed_expression_;
		if ( last_constructed_expression_ == nullptr ) {
			utility_exit_with_message( "Error constructing expression objects for ArithmeticASTFunction: argument #" +
				utility::to_string( node.function()->nargs() ) + " to function: " + node.function()->name() + " resulted in NULL pointer when parsed" );
		}
	}
	/// Polymorphic dispatch to allow derived classes to intercept the handling of this function.
	last_constructed_expression_ = handle_function_expression( node.function(), args );

	if ( last_constructed_expression_ == nullptr ) {
		/// ERROR.
		utility_exit_with_message( "ExpressionCreator::visit( ArithmeticASTFunction ) called with unrecognized function: " +
			node.function()->name() + ".\n" +
			"Drived classes also failed to recognize this function.");
	}
}

void
ExpressionCreator::visit( ArithmeticASTTerm const & node )
{
	utility::vector1< ExpressionCOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		last_constructed_expression_.reset(); // in case the child doesn't zero this out?
		(*iter)->visit( *this );
		if ( last_constructed_expression_ ) {
			expressions.push_back( last_constructed_expression_ );
			semi_constructed_expression_ = last_constructed_expression_;
		}
	}
	if ( expressions.size() == 1 ) {
		last_constructed_expression_ = expressions[ 1 ];
		semi_constructed_expression_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_expression_ = expressions[ 2 ];
		semi_constructed_expression_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTTerm.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}
}

void
ExpressionCreator::visit( ArithmeticASTFactor const & node )
{
	node.child()->visit( *this );
}

void
ExpressionCreator::visit( ArithmeticASTValue const & node )
{
	if ( node.is_literal() ) {
		last_constructed_expression_ = ExpressionCOP( ExpressionOP( new LiteralExpression( node.literal_value() ) ) );
	} else {
		last_constructed_expression_ = handle_variable_expression( node );
	}
}

void
ExpressionCreator::visit( ArithmeticASTRestTerm const & node )
{
	if ( node.children_begin() == node.children_end() ) {
		last_constructed_expression_.reset();
		semi_constructed_expression_.reset();
	} else {
		ExpressionCOP parents_semi_constructed_expression = semi_constructed_expression_;
		ExpressionCOP my_completed_expression;
		semi_constructed_expression_.reset();

		utility::vector1< ExpressionCOP > expressions;
		expressions.reserve( 2 );

		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			last_constructed_expression_.reset(); // in case the child doesn't zero this out?
			(*iter)->visit( *this );
			if ( last_constructed_expression_ ) {
				expressions.push_back( last_constructed_expression_ );
				if ( iter == node.children_begin() ) {

					if ( node.rest_term_token() == MULTIPLY_SYMBOL ) {
						semi_constructed_expression_ = ExpressionCOP( ExpressionOP( new MultiplyExpression( parents_semi_constructed_expression, last_constructed_expression_ ) ) );
					} else if ( node.rest_term_token() == DIVIDE_SYMBOL ) {
						semi_constructed_expression_ = ExpressionCOP( ExpressionOP( new DivideExpression( parents_semi_constructed_expression, last_constructed_expression_ ) ) );
					} else {
						utility_exit_with_message( "Error visiting ArithmeticASTRestExpression: expected MULTIPLY_SYMBOL or DIVIDE_SYMBOL; got " + token_type_name( node.rest_term_token() ) );
					}

					my_completed_expression = semi_constructed_expression_;
				}
			}
		}
		if ( expressions.size() == 1 ) {
			last_constructed_expression_ = my_completed_expression;
		} else if ( expressions.size() == 2 ) {
			// last_constructed_expression_ is already current and was set by child
			assert( last_constructed_expression_ );
			semi_constructed_expression_.reset();
		} else {
			utility_exit_with_message( "Error in visiting children of ArithmeticASTRestTerm: too many children ("
				+ utility::to_string( expressions.size() ) + ")" );
		}
	}
}

void
ExpressionCreator::visit( ArithmeticASTRestExpression const & node )
{
	if ( node.children_begin() == node.children_end() ) {
		last_constructed_expression_.reset();
		semi_constructed_expression_.reset();
	} else {
		ExpressionCOP parents_semi_constructed_expression = semi_constructed_expression_;
		ExpressionCOP my_completed_expression;
		semi_constructed_expression_.reset();

		utility::vector1< ExpressionCOP > expressions;
		expressions.reserve( 2 );

		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			last_constructed_expression_.reset(); // in case the child doesn't zero this out?
			(*iter)->visit( *this );
			if ( last_constructed_expression_ ) {
				expressions.push_back( last_constructed_expression_ );
				if ( iter == node.children_begin() ) {
					if ( node.rest_expression_token() == PLUS_SYMBOL ) {
						semi_constructed_expression_ = ExpressionCOP( ExpressionOP( new AddExpression( parents_semi_constructed_expression, last_constructed_expression_ ) ) );
					} else if ( node.rest_expression_token() == SUBTRACT_SYMBOL ) {
						semi_constructed_expression_ = ExpressionCOP( ExpressionOP( new SubtractExpression( parents_semi_constructed_expression, last_constructed_expression_ ) ) );
					} else {
						utility_exit_with_message( "Error visiting ArithmeticASTRestTerm: expected PLUS_SYMBOL or SUBTRACT_SYMBOL; got " + token_type_name( node.rest_expression_token() ) );
					}
					my_completed_expression = semi_constructed_expression_;
				}
			}
		}
		if ( expressions.size() == 1 ) {
			last_constructed_expression_ = my_completed_expression;
		} else if ( expressions.size() == 2 ) {
			// last_constructed_expression_ is already current and was set by child
			assert( last_constructed_expression_ );
			semi_constructed_expression_.reset();
		} else {
			utility_exit_with_message( "Error in visiting children of ArithmeticASTRestExpression: too many children ("
				+ utility::to_string( expressions.size() ) + ")" );
		}
	}
}


void
ExpressionCreator::visit( ArithmeticASTNode const & )
{
	utility_exit_with_message( "Dispatch to unrecognized ArithmeticASTNode type in ExpressionCreator!" );
}

ExpressionCOP
ExpressionCreator::create_expression_tree( ArithmeticASTExpression const & expr )
{
	last_constructed_expression_.reset();
	semi_constructed_expression_.reset();
	visit( expr );
	return last_constructed_expression_;
}

/// @brief Factory method to be implemented by derived classes
/// which may wish to handle variable expressions differently.
ExpressionCOP
ExpressionCreator::handle_variable_expression( ArithmeticASTValue const & )
{
	utility_exit_with_message( "ExpressionCreator::handle_variable_expression called.\nThis  method should be overridden by the derived class.");
	return nullptr;
}

ExpressionCOP
ExpressionCreator::handle_function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
)
{
	Size count_args = args.size();
	std::string fname = function->name();
	if ( fname == "max" ) {
		if ( count_args != 2 ) {
			utility_exit_with_message( "Error constructing MaxExpression; Did not create 2 arguments, created "
				+ utility::to_string( count_args ) + " argument" + (count_args == 1 ? "." : "s.") );
		}
		return ExpressionCOP( ExpressionOP( new MaxExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "min" ) {
		if ( count_args != 2 ) {
			utility_exit_with_message( "Error constructing MinExpression; Did not create 2 arguments, created "
				+ utility::to_string( count_args ) +" argument"+ (count_args == 1 ? "." : "s.") );
		}
		return ExpressionCOP( ExpressionOP( new MinExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "sqrt" ) {
		if ( count_args != 1 ) {
			utility_exit_with_message( "Error constructing SquarerootExpression; Did not create 1 arguments, created "
				+ utility::to_string( count_args ) +" arguments." );
		}
		return ExpressionCOP( ExpressionOP( new SquarerootExpression( args[ 1 ] ) ) );
	}

	return nullptr;
}

SimpleExpressionCreator::SimpleExpressionCreator() = default;

SimpleExpressionCreator::SimpleExpressionCreator(
	std::list< std::string > const & varnames
)
{
	for ( auto const & varname : varnames ) {
		add_variable( varname );
	}
}


void
SimpleExpressionCreator::add_variable( std::string const & varname )
{
	if ( variables_.find( varname ) != variables_.end() ) {
		utility_exit_with_message( "Error adding variable: '" + varname + "'; already present in variables_ map." );
	}
	variables_[ varname ] = VariableExpressionOP( new VariableExpression( varname ) );
}


ExpressionCOP
SimpleExpressionCreator::handle_variable_expression( ArithmeticASTValue const & node )
{
	if ( node.is_literal() ) {
		utility_exit_with_message( "Error in SimpleExpressionCreator::handle_variable_expression; non-variable (literal) node given!" +
			utility::to_string( node.literal_value() ));
	}
	if ( variables_.find( node.variable_name() ) == variables_.end() ) {
		utility_exit_with_message( "Error in SimpleExpressionCreator::handle_variable_expression; variable '"
			+ node.variable_name() + "' not found in variables_ map" );
	}
	return variables_.find( node.variable_name() )->second;
}

VariableExpressionOP
SimpleExpressionCreator::get_variable( std::string const & varname )
{
	if ( variables_.find( varname ) == variables_.end() ) {
		utility_exit_with_message( "Error in SimpleExpressionCreator::get_variable; variable '"
			+ varname + "' not found in variables_ map" );
	}
	return variables_.find( varname )->second;

}

std::map< std::string, VariableExpressionOP >
SimpleExpressionCreator::variables() const
{
	return variables_;
}


BooleanExpressionCreator::BooleanExpressionCreator() : SimpleExpressionCreator() {}
BooleanExpressionCreator::BooleanExpressionCreator( std::list< std::string > const & varnames )
:
	SimpleExpressionCreator( varnames )
{}

ExpressionCOP
BooleanExpressionCreator::handle_function_expression(
	FunctionTokenCOP function,
	utility::vector1< ExpressionCOP > const & args
)
{
	ExpressionCOP parent_result = grandparent::handle_function_expression( function, args );
	if ( parent_result ) return parent_result;

	std::string const fname = function->name();
	if ( fname == "EQUALS" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new EqualsExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "GT" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new GT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "GTE" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new GTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "LT" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new LT_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "LTE" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new LTE_Expression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname ==  "AND" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new AndExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "OR" ) {
		assert( args.size() == 2 );
		return ExpressionCOP( ExpressionOP( new OrExpression( args[ 1 ], args[ 2 ] ) ) );
	} else if ( fname == "NOT" ) {
		assert( args.size() == 1 );
		return ExpressionCOP( ExpressionOP( new NotExpression( args[ 1 ] ) ) );
	} else {
		utility_exit_with_message( "Unrecognized function name in BooleanExpressionCreator: '" + fname + "'" );
		return nullptr;
	}
}


LiteralExpression::LiteralExpression() : value_( 0.0 )
{}

LiteralExpression::LiteralExpression( numeric::Real value ) : value_( value )
{}

void LiteralExpression::set_value( numeric::Real value )
{
	value_ = value;
}

numeric::Real
LiteralExpression::operator() () const
{
	return value_;
}


ExpressionCOP
LiteralExpression::differentiate( std::string const & ) const
{
	return nullptr;
}

std::list< std::string >
LiteralExpression::active_variables() const
{
	std::list< std::string > empty;
	return empty;
}

VariableExpression::VariableExpression( std::string  name )
:
	name_(std::move( name )),
	value_( 0 )
{}

VariableExpression::VariableExpression( std::string  name, numeric::Real value )
:
	name_(std::move( name )),
	value_( value )
{}


numeric::Real
VariableExpression::operator() () const
{
	return value_;
}


void VariableExpression::set_value( numeric::Real value )
{
	value_ = value;
}


std::string
VariableExpression::name() const
{
	return name_;
}

ExpressionCOP
VariableExpression::differentiate( std::string const & varname ) const
{
	if ( name_ == varname ) {
		return ExpressionCOP( ExpressionOP( new LiteralExpression( 1.0 ) ) );
	}
	return nullptr;
}

std::list< std::string >
VariableExpression::active_variables() const
{
	std::list< std::string > myname;
	myname.push_back( name_ );
	return myname;
}


UnaryExpression::UnaryExpression() : ex_( /* 0 */ ) {}
UnaryExpression::~UnaryExpression() = default;

UnaryExpression::UnaryExpression( ExpressionCOP ex ) : ex_(std::move( ex )) {}

void
UnaryExpression::set_expression( ExpressionCOP ex ) { ex_ = ex; }

ExpressionCOP
UnaryExpression::ex() const {
	if ( ex_ == nullptr ) {
		utility_exit_with_message( "Bad expression evaluation; protocols::optimize_weights::UnaryExpression::ex_ is null" );
	}
	return ex_;
}

std::list< std::string >
UnaryExpression::active_variables() const
{
	return ex_->active_variables();
}


BinaryExpression::BinaryExpression() : e1_( /* 0 */ ), e2_( nullptr ) {}
BinaryExpression::~BinaryExpression() = default;
BinaryExpression::BinaryExpression( ExpressionCOP e1, ExpressionCOP e2 ) : e1_(std::move( e1 )), e2_(std::move( e2 )) {}

void BinaryExpression::set_first_expression( ExpressionCOP e1 ) { e1_ = e1; }
void BinaryExpression::set_second_expression( ExpressionCOP e2 ) { e2_ = e2; }

ExpressionCOP
BinaryExpression::e1() const
{
	if ( e1_ == nullptr ) {
		utility_exit_with_message( "Bad expression evaluation; protocols::optimize_weights::BinaryExpression::e1_ is null" );
	}
	return e1_;
}

ExpressionCOP
BinaryExpression::e2() const
{
	if ( e2_ == nullptr ) {
		utility_exit_with_message( "Bad expression evaluation; protocols::optimize_weights::BinaryExpression::e2_ is null" );
	}
	return e2_;
}

std::list< std::string >
BinaryExpression::active_variables() const
{
	std::list< std::string > childrens_vars;
	std::list< std::string > e1_vars = e1_->active_variables();
	std::list< std::string > e2_vars = e2_->active_variables();
	childrens_vars.splice( childrens_vars.end(), e1_vars );
	childrens_vars.splice( childrens_vars.end(), e2_vars );
	return childrens_vars;
}

SquarerootExpression::SquarerootExpression() : UnaryExpression() {}
SquarerootExpression::SquarerootExpression( ExpressionCOP ex ) : UnaryExpression( ex ) {}

numeric::Real
SquarerootExpression::operator() () const
{
	return std::sqrt( (*ex())() );
}

/// @details Badly behaved expression around 0
ExpressionCOP
SquarerootExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP dex_dvar = ex()->differentiate( varname );
	if ( dex_dvar ) {
		LiteralExpressionOP onehalf( new LiteralExpression( 0.5 ) );
		//LiteralExpressionOP one = new LiteralExpression( 1.0 );
		SquarerootExpressionOP sqrt( new SquarerootExpression( ex() ) );
		DivideExpressionOP invsqrt( new DivideExpression( dex_dvar, sqrt ) );
		MultiplyExpressionOP derivative( new MultiplyExpression( onehalf, invsqrt ) );
		return derivative;
	}
	return nullptr;
}

AbsoluteValueExpression::AbsoluteValueExpression() : UnaryExpression() {}
AbsoluteValueExpression::AbsoluteValueExpression( ExpressionCOP ex ) : UnaryExpression( ex ) {}

numeric::Real
AbsoluteValueExpression::operator() () const
{
	return std::abs( (*ex())() );
}

/// @details Badly behaved expression around 0
ExpressionCOP
AbsoluteValueExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP dex_dvar = ex()->differentiate( varname );

	if ( ! dex_dvar ) return nullptr;

	if (  (*ex())() > 0 ) {
		return dex_dvar;
	} else {
		LiteralExpressionOP negone( new LiteralExpression( -1 ) );
		return ExpressionCOP( ExpressionOP( new MultiplyExpression( negone, dex_dvar ) ) );
	}
}

AddExpression::AddExpression() : BinaryExpression() {}
AddExpression::AddExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides of the expression; return the sum.
numeric::Real
AddExpression::operator() () const
{
	return (*e1())() + (*e2())();
}

ExpressionCOP
AddExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	} else if ( de1 && de2 ) {
		return ExpressionCOP( ExpressionOP( new AddExpression( de1, de2 ) ) );
	} else if ( de1 ) {
		return de1;
	} else {
		return de2;
	}

}

SubtractExpression::SubtractExpression() : BinaryExpression() {}
SubtractExpression::SubtractExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides of the expression; return the sum.
numeric::Real
SubtractExpression::operator() () const
{
	return (*e1())() - (*e2())();
}

ExpressionCOP
SubtractExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	} else if ( de1 && de2 ) {
		return ExpressionCOP( ExpressionOP( new SubtractExpression( de1, de2 ) ) );
	} else if ( de1 ) {
		return de1;
	} else {
		LiteralExpressionCOP negone( LiteralExpressionOP( new LiteralExpression( -1 ) ) );
		return ExpressionCOP( ExpressionOP( new MultiplyExpression( negone, de2 ) ) );
	}

}

MultiplyExpression::MultiplyExpression() : BinaryExpression() {}
MultiplyExpression::MultiplyExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides; return the product.
numeric::Real
MultiplyExpression::operator() () const
{
	return (*e1())() * (*e2())();
}

ExpressionCOP
MultiplyExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	} else if ( de1 && de2 ) {
		ExpressionCOP a( ExpressionOP( new MultiplyExpression( de1, e2() ) ) );
		ExpressionCOP b( ExpressionOP( new MultiplyExpression( e1(), de2 ) ) );
		return ExpressionCOP( ExpressionOP( new AddExpression( a, b ) ) );
	} else if ( de1 ) {
		return ExpressionCOP( ExpressionOP( new MultiplyExpression( de1, e2() ) ) );
	} else {
		return ExpressionCOP( ExpressionOP( new MultiplyExpression( de2, e1() ) ) );
	}

}

DivideExpression::DivideExpression() : BinaryExpression() {}
DivideExpression::DivideExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides; return the product.
numeric::Real
DivideExpression::operator() () const
{
	return (*e1())() / (*e2())();
}

ExpressionCOP
DivideExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	} else if ( de1 && de2 ) {
		MultiplyExpressionOP num1( new MultiplyExpression( de1, e2() ) );
		MultiplyExpressionOP num2( new MultiplyExpression( de2, e1() ) );
		SubtractExpressionOP diff( new SubtractExpression( num1, num2 ) );
		MultiplyExpressionOP sqr( new MultiplyExpression( e2(), e2() ) );
		return ExpressionCOP( ExpressionOP( new DivideExpression( diff, sqr ) ) );
	} else if ( de1 ) {
		return ExpressionCOP( ExpressionOP( new DivideExpression( de1, e2() ) ) );
	} else {
		LiteralExpressionOP negone( new LiteralExpression( -1.0 ) );
		MultiplyExpressionOP de2_e1( new MultiplyExpression( de2, e1() ) );
		MultiplyExpressionOP num( new MultiplyExpression( negone, de2_e1 ) );
		MultiplyExpressionOP sqr( new MultiplyExpression( e2(), e2() ) );
		return ExpressionCOP( ExpressionOP( new DivideExpression( num, sqr ) ) );
	}

}

MaxExpression::MaxExpression() : BinaryExpression() {}
MaxExpression::MaxExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides; return the product.
numeric::Real
MaxExpression::operator() () const
{
	return std::max((*e1())(), (*e2())());
}

ExpressionCOP
MaxExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	}

	if ( de1 == nullptr ) de1 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	if ( de2 == nullptr ) de2 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );

	return ExpressionCOP( ExpressionOP( new MetaMaxExpression( e1(), e2(), de1, de2 ) ) );

}

std::list< std::string >
MaxExpression::active_variables() const
{
	return (*e1())() > (*e2())() ? e1()->active_variables() : e2()->active_variables();
}


MinExpression::MinExpression() : BinaryExpression() {}
MinExpression::MinExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}

/// @details get expression pointers; evaluate left and right hand sides; return the product.
numeric::Real
MinExpression::operator() () const
{
	return std::min((*e1())(), (*e2())());
}

ExpressionCOP
MinExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP de1 = e1()->differentiate( varname );
	ExpressionCOP de2 = e2()->differentiate( varname );

	if ( ! de1 && ! de2 ) {
		return nullptr;
	}

	if ( de1 == nullptr ) de1 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	if ( de2 == nullptr ) de2 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );

	return ExpressionCOP( ExpressionOP( new MetaMinExpression( e1(), e2(), de1, de2 ) ) );

}

std::list< std::string >
MinExpression::active_variables() const
{
	return (*e1())() < (*e2())() ? e1()->active_variables() : e2()->active_variables();
}


/// @brief Evaluates ee1 when e1 is larger than e2; evaluates ee2 otherwise.
MetaMaxExpression::MetaMaxExpression(
	ExpressionCOP e1,
	ExpressionCOP e2,
	ExpressionCOP ee1,
	ExpressionCOP ee2
) :
	e1_(std::move( e1 )),
	e2_(std::move( e2 )),
	ee1_(std::move( ee1 )),
	ee2_(std::move( ee2 ))
{}

numeric::Real
MetaMaxExpression::operator() () const
{
	return (*e1_)() > (*e2_)() ? (*ee1_)() : (*ee2_ )();
}

ExpressionCOP
MetaMaxExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP dee1 = ee1_->differentiate( varname );
	ExpressionCOP dee2 = ee2_->differentiate( varname );

	if ( ! dee1 && ! dee2 ) return nullptr;

	if ( ! dee1 ) dee1 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	if ( ! dee2 ) dee2 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	return ExpressionCOP( ExpressionOP( new MetaMaxExpression( e1_, e2_, dee1, dee2 ) ) );

}

std::list< std::string >
MetaMaxExpression::active_variables() const
{
	return (*e1_)() > (*e2_)() ? ee1_->active_variables() : ee2_->active_variables();
}


MetaMinExpression::MetaMinExpression(
	ExpressionCOP e1,
	ExpressionCOP e2,
	ExpressionCOP ee1,
	ExpressionCOP ee2
) :
	e1_(std::move( e1 )),
	e2_(std::move( e2 )),
	ee1_(std::move( ee1 )),
	ee2_(std::move( ee2 ))
{}

numeric::Real
MetaMinExpression::operator() () const
{
	return (*e1_)() < (*e2_)() ? (*ee1_)() : (*ee2_ )();
}

ExpressionCOP
MetaMinExpression::differentiate( std::string const & varname ) const
{
	ExpressionCOP dee1 = ee1_->differentiate( varname );
	ExpressionCOP dee2 = ee2_->differentiate( varname );

	if ( ! dee1 && ! dee2 ) return nullptr;

	if ( ! dee1 ) dee1 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	if ( ! dee2 ) dee2 = ExpressionCOP( ExpressionOP( new LiteralExpression( 0.0 ) ) );
	return ExpressionCOP( ExpressionOP( new MetaMinExpression( e1_, e2_, dee1, dee2 ) ) );

}

std::list< std::string >
MetaMinExpression::active_variables() const
{
	return (*e1_)() < (*e2_)() ? ee1_->active_variables() : ee2_->active_variables();
}


EqualsExpression::EqualsExpression() : BinaryExpression() {}
EqualsExpression::EqualsExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
EqualsExpression::~EqualsExpression() = default;


numeric::Real
EqualsExpression::operator() () const
{
	return (*e1())() == (*e2())();
}


ExpressionCOP
EqualsExpression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


GT_Expression::GT_Expression() : BinaryExpression() {}
GT_Expression::GT_Expression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
GT_Expression::~GT_Expression() = default;


numeric::Real
GT_Expression::operator() () const
{
	return (*e1())() > (*e2())();
}


ExpressionCOP
GT_Expression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


GTE_Expression::GTE_Expression() : BinaryExpression() {}
GTE_Expression::GTE_Expression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
GTE_Expression::~GTE_Expression() = default;


numeric::Real
GTE_Expression::operator() () const
{
	return (*e1())() >= (*e2())();
}


ExpressionCOP
GTE_Expression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


LT_Expression::LT_Expression() : BinaryExpression() {}
LT_Expression::LT_Expression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
LT_Expression::~LT_Expression() = default;


numeric::Real
LT_Expression::operator() () const
{
	return (*e1())() < (*e2())();
}


ExpressionCOP
LT_Expression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


LTE_Expression::LTE_Expression() : BinaryExpression() {}
LTE_Expression::LTE_Expression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
LTE_Expression::~LTE_Expression() = default;


numeric::Real
LTE_Expression::operator() () const
{
	return (*e1())() <= (*e2())();
}


ExpressionCOP
LTE_Expression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


/// 2. Boolean Logic Operators
AndExpression::AndExpression() : BinaryExpression() {}
AndExpression::AndExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
AndExpression::~AndExpression() = default;


numeric::Real
AndExpression::operator() () const
{
	return (*e1())() != 0.0 && (*e2())() != 0.0;
}


ExpressionCOP
AndExpression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


OrExpression::OrExpression() : BinaryExpression() {}
OrExpression::OrExpression( ExpressionCOP e1, ExpressionCOP e2 ) : BinaryExpression( e1, e2 ) {}
OrExpression::~OrExpression() = default;

numeric::Real
OrExpression::operator() () const
{
	return (*e1())() != 0.0 || (*e2())() != 0.0;
}


ExpressionCOP
OrExpression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


NotExpression::NotExpression() : UnaryExpression() {}
NotExpression::NotExpression( ExpressionCOP e ) : UnaryExpression( e ) {}
NotExpression::~NotExpression() = default;

numeric::Real
NotExpression::operator() () const
{
	return (*ex())() == 0.0 ? 1.0 : 0.0;
}

ExpressionCOP
NotExpression::differentiate( std::string const & ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}


ITEExpression::ITEExpression(
	ExpressionCOP condition,
	ExpressionCOP then_expression,
	ExpressionCOP else_expression
) :
	condition_(std::move( condition )),
	then_expression_(std::move( then_expression )),
	else_expression_(std::move( else_expression ))
{}

numeric::Real
ITEExpression::operator() () const
{
	numeric::Real condition_val = (*condition_)();
	//std::cout << "ITE(); con=" << condition_val << "ret=";
	if ( condition_val == 0.0 ) { // 0.0 is false, everything else is true.
		//std::cout << (*else_expression_)() << std::endl;
		return (*else_expression_)();
	} else {
		//std::cout << (*then_expression_)() << std::endl;
		return (*then_expression_)();
	}
}

ExpressionCOP
ITEExpression::differentiate( std::string const &  ) const
{
	/// Boolean expressions cannot be differentiated
	return nullptr;
}

ExpressionCOP
ITEExpression::condition() const {
	return condition_;
}

ExpressionCOP
ITEExpression::then_expression() const {
	return then_expression_;
}

ExpressionCOP
ITEExpression::else_expression() const {
	return else_expression_;
}


std::list< std::string >
ITEExpression::active_variables() const
{
	numeric::Real condition_val = (*condition_)();
	if ( condition_val == 0.0 ) { // 0.0 is false, everything else is true.
		return else_expression_->active_variables();
	} else {
		return then_expression_->active_variables();
	}
}


ExpressionCOP
parse_string_to_expression( std::string const & input_string )
{
	ArithmeticScanner as;
	TokenSetOP tokens = as.scan( input_string );
	ArithmeticASTExpression expr_ast;
	expr_ast.parse( *tokens );
	ExpressionCreator ec;
	return ec.create_expression_tree( expr_ast );
}

ExpressionCOP
parse_string_to_boolean_expression( std::string const & input_string )
{
	BooleanExpressionScanner bas;
	TokenSetOP tokens = bas.scan( input_string );
	ArithmeticASTExpression expr_ast;
	expr_ast.parse( *tokens );
	BooleanExpressionCreator ec;
	return ec.create_expression_tree( expr_ast );
}

}
}

