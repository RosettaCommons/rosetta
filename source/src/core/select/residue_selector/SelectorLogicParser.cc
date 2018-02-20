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

// Unit headers
#include <core/select/residue_selector/SelectorLogicParser.hh>

// Package headers
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <iostream>

namespace core {
namespace select {
namespace residue_selector {

// Expression visitor
SelectorLogicParser::SelectorLogicParser() = default;

SelectorLogicParser::~SelectorLogicParser() = default;

ResidueSelectorOP
SelectorLogicParser:: parse_string_to_residue_selector(
	basic::datacache::DataMap const & dm,
	std::string const & selector_logic_string
)
{
	if ( ! dm.has_type( "ResidueSelector" ) ) {
		std::ostringstream oss;
		oss << "There are no residue selectors in the DataMap; have any been declared?\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	// create scanner w/o adding any functions, e.g. sqrt
	numeric::expression_parser::ArithmeticScanner as( true );
	as.treat_AND_and_OR_as_operators( true );

	for ( auto const & name_selector_pair : dm.category_map( "ResidueSelector" ) ) {
		ResidueSelectorOP selector = utility::pointer::dynamic_pointer_cast< ResidueSelector >
			(  name_selector_pair.second );
		if ( selector ) {
			as.add_variable( name_selector_pair.first );
			add_named_selector( name_selector_pair.first, selector );
		}
	}
	return finish_parse_string_to_residue_selector( selector_logic_string, as );
}


ResidueSelectorOP
SelectorLogicParser::parse_string_to_residue_selector(
	std::string const & input_string,
	std::map< std::string, ResidueSelectorOP > const & selectors
)
{
	// create scanner w/o adding any functions, e.g. sqrt
	numeric::expression_parser::ArithmeticScanner as( true );
	as.treat_AND_and_OR_as_operators( true );

	for ( auto const & name_selector_pair : selectors ) {
		as.add_variable( name_selector_pair.first );
		add_named_selector( name_selector_pair.first, name_selector_pair.second );
	}
	return finish_parse_string_to_residue_selector( input_string, as );

}

void SelectorLogicParser::add_named_selector(
	std::string const & name,
	ResidueSelectorOP selector
)
{
	// TO DO: verify no repeats
	variables_[ name ] = selector;
}


void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTExpression const & node)
{
	utility::vector1< ResidueSelectorOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		if ( last_constructed_selector_ ) {
			expressions.push_back( last_constructed_selector_ );
			semi_constructed_selector_ = last_constructed_selector_;
		}
	}

	if ( expressions.size() == 1 ) {
		last_constructed_selector_ = expressions[ 1 ];
		semi_constructed_selector_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_selector_ = expressions[ 2 ];
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTExpression.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}

}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTOrClause const & node)
{
	utility::vector1< ResidueSelectorOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		if ( last_constructed_selector_ ) {
			expressions.push_back( last_constructed_selector_ );
			semi_constructed_selector_ = last_constructed_selector_;
		}
	}

	if ( expressions.size() == 1 ) {
		last_constructed_selector_ = expressions[ 1 ];
		semi_constructed_selector_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_selector_ = expressions[ 2 ];
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTOrClause.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}

}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTAndClause const & node)
{
	utility::vector1< ResidueSelectorOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		(*iter)->visit( *this );
		if ( last_constructed_selector_ ) {
			expressions.push_back( last_constructed_selector_ );
			semi_constructed_selector_ = last_constructed_selector_;
		}
	}

	if ( expressions.size() == 1 ) {
		last_constructed_selector_ = expressions[ 1 ];
		semi_constructed_selector_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_selector_ = expressions[ 2 ];
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTAndClause.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}

}


void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTFunction const & node)
{
	utility::vector1< ResidueSelectorOP > args( node.function()->nargs(), nullptr );
	Size count_args = 0;
	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		++count_args;
		(*iter)->visit( *this );
		args[ count_args ] = last_constructed_selector_;
		if ( last_constructed_selector_ == nullptr ) {
			utility_exit_with_message( "Error constructing expression objects for ArithmeticASTFunction: argument #" +
				utility::to_string( node.function()->nargs() ) + " to function: " + node.function()->name() + " resulted in NULL pointer when parsed" );
		}
	}
	/// Polymorphic dispatch to allow derived classes to intercept the handling of this function.
	last_constructed_selector_ = handle_function_expression( node.function(), args );

	if ( last_constructed_selector_ == nullptr ) {
		/// ERROR.
		utility_exit_with_message( "SelectorLogicParser::visit( ArithmeticASTFunction ) called with unrecognized function: " +
			node.function()->name() + ".\n" +
			"Drived classes also failed to recognize this function.");
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTTerm const & node )
{
	utility::vector1< ResidueSelectorOP > expressions;
	expressions.reserve( 2 );

	for ( auto iter = node.children_begin(),
			iter_end = node.children_end(); iter != iter_end; ++iter ) {
		last_constructed_selector_.reset(); // in case the child doesn't zero this out?
		(*iter)->visit( *this );
		if ( last_constructed_selector_ ) {
			expressions.push_back( last_constructed_selector_ );
			semi_constructed_selector_ = last_constructed_selector_;
		}
	}
	if ( expressions.size() == 1 ) {
		last_constructed_selector_ = expressions[ 1 ];
		semi_constructed_selector_.reset();
	} else if ( expressions.size() == 2 ) {
		last_constructed_selector_ = expressions[ 2 ];
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Visited too many children of ArithmeticASTTerm.  Cannot proceed! #children= " +
			utility::to_string( expressions.size() ) );
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTFactor const & node )
{
	if ( node.not_token() == numeric::expression_parser::NOT_SYMBOL ) {
		node.child()->visit( *this );
		ResidueSelectorOP not_selector( new NotResidueSelector( last_constructed_selector_ ) );
		last_constructed_selector_ = not_selector;
		semi_constructed_selector_.reset();
	} else {
		node.child()->visit( *this );
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTValue const & node )
{
	if ( node.is_literal() ) {
		utility_exit_with_message( "Cannot handle literals (e.g. \"4\") when creating logical ResidueSelectors from strings" );
	} else {
		last_constructed_selector_ = handle_variable_expression( node );
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTRestTerm const & node )
{
	if ( node.children_begin() == node.children_end() ) {
		last_constructed_selector_.reset();
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Error visiting ArithmeticASTRestExpression: Symbol " + token_type_name( node.rest_term_token() ) + " cannot be processed when constructing a ResidueSelector" );
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTRestAndClause const & node )
{
	if ( node.children_begin() == node.children_end() ) {
		last_constructed_selector_.reset();
		semi_constructed_selector_.reset();
	} else {
		utility_exit_with_message( "Error visiting ArithmeticASTRestAndClause: Symbol " + token_type_name( node.rest_and_clause_token() ) + " cannot be processed when constructing a ResidueSelector" );
	}
}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTRestOrClause const & node )
{

	if ( node.children_begin() == node.children_end() ) {
		last_constructed_selector_.reset();
		semi_constructed_selector_.reset();
	} else {
		ResidueSelectorOP parents_semi_constructed_expression = semi_constructed_selector_;
		ResidueSelectorOP my_completed_expression;
		semi_constructed_selector_.reset();

		utility::vector1< ResidueSelectorOP > expressions;
		expressions.reserve( 2 );

		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			last_constructed_selector_.reset(); // in case the child doesn't zero this out?
			(*iter)->visit( *this );
			if ( last_constructed_selector_ ) {
				expressions.push_back( last_constructed_selector_ );
				if ( iter == node.children_begin() ) {
					if ( node.rest_or_clause_token() == numeric::expression_parser::AND_SYMBOL ) {
						semi_constructed_selector_.reset( new AndResidueSelector( parents_semi_constructed_expression, last_constructed_selector_ ) );
					} else {
						utility_exit_with_message( "Error visiting ArithmeticASTRestOrClause: expected AND_SYMBOL; got " + token_type_name( node.rest_or_clause_token() ) );
					}
					my_completed_expression = semi_constructed_selector_;
				}
			}
		}
		if ( expressions.size() == 1 ) {
			last_constructed_selector_ = my_completed_expression;
		} else if ( expressions.size() == 2 ) {
			// last_constructed_selector_ is already current and was set by child
			assert( last_constructed_selector_ );
			semi_constructed_selector_.reset();
		} else {
			utility_exit_with_message( "Error in visiting children of ArithmeticASTRestOrClause: too many children ("
				+ utility::to_string( expressions.size() ) + ")" );
		}
	}

}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTRestExpression const & node )
{

	if ( node.children_begin() == node.children_end() ) {
		last_constructed_selector_.reset();
		semi_constructed_selector_.reset();
	} else {
		ResidueSelectorOP parents_semi_constructed_expression = semi_constructed_selector_;
		ResidueSelectorOP my_completed_expression;
		semi_constructed_selector_.reset();

		utility::vector1< ResidueSelectorOP > expressions;
		expressions.reserve( 2 );

		for ( auto iter = node.children_begin(),
				iter_end = node.children_end(); iter != iter_end; ++iter ) {
			last_constructed_selector_.reset(); // in case the child doesn't zero this out?
			(*iter)->visit( *this );
			if ( last_constructed_selector_ ) {
				expressions.push_back( last_constructed_selector_ );
				if ( iter == node.children_begin() ) {
					if ( node.rest_expression_token() == numeric::expression_parser::OR_SYMBOL ) {
						semi_constructed_selector_.reset( new OrResidueSelector( parents_semi_constructed_expression, last_constructed_selector_ ) );
					} else {
						utility_exit_with_message( "Error visiting ArithmeticASTRestExpression: expected OR_SYMBOL; got " + token_type_name( node.rest_expression_token() ) );
					}
					my_completed_expression = semi_constructed_selector_;
				}
			}
		}
		if ( expressions.size() == 1 ) {
			last_constructed_selector_ = my_completed_expression;
		} else if ( expressions.size() == 2 ) {
			// last_constructed_selector_ is already current and was set by child
			assert( last_constructed_selector_ );
			semi_constructed_selector_.reset();
		} else {
			utility_exit_with_message( "Error in visiting children of ArithmeticASTRestExpression: too many children ("
				+ utility::to_string( expressions.size() ) + ")" );
		}
	}

}

void
SelectorLogicParser::visit( numeric::expression_parser::ArithmeticASTNode const & )
{
	utility_exit_with_message( "Dispatch to unrecognized ArithmeticASTNode type in SelectorLogicParser!" );
}

ResidueSelectorOP
SelectorLogicParser::create_selector(
	numeric::expression_parser::ArithmeticASTExpression const & expr
)
{
	last_constructed_selector_.reset();
	semi_constructed_selector_.reset();
	visit( expr );
	return last_constructed_selector_;
}

/// @brief Factory method to be implemented by derived classes
/// which may wish to handle variable expressions differently.
ResidueSelectorOP
SelectorLogicParser::handle_variable_expression( numeric::expression_parser::ArithmeticASTValue const & node )
{
	auto iter = variables_.find( node.variable_name() );
	if ( iter == variables_.end() ) {
		utility_exit_with_message( "Error requesting undefined ResidueSelector with name \"" +
			node.variable_name() + "\"; this error should not have been reached, as any un-detected" +
			" variables should have been caught at the tokenization stage?" );
	}
	return iter->second;
}

ResidueSelectorOP
SelectorLogicParser::handle_function_expression(
	numeric::expression_parser::FunctionTokenCOP,
	utility::vector1< ResidueSelectorOP > const &
)
{
	utility_exit_with_message( "Currently no support for defining ResidueSelector functions" );
}


ResidueSelectorOP
SelectorLogicParser::finish_parse_string_to_residue_selector(
	std::string const & input_string,
	numeric::expression_parser::ArithmeticScanner & as
)
{
	using namespace numeric::expression_parser;

	TokenSetOP tokens;
	try {
		tokens = as.scan( input_string );
	} catch ( utility::excn::Exception & e ) {
		std::ostringstream oss;
		oss << "Failed to tokenize string logically combining ResidueSelectors.\n";
		oss << "\"" << input_string << "\"\n";
		oss << "Allowed names for ResidueSelectors are:\n";
		for ( auto const & name_selector_pair : variables_ ) {
			oss << "   " << name_selector_pair.first << "\n";
		}
		oss << "Error message from Scanner:\n";
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	ArithmeticASTExpression expr_ast;
	try {
		expr_ast.parse( *tokens );
	} catch ( utility::excn::Exception & e ) {
		std::ostringstream oss;
		oss << "Failed to parse the string logically combining ResidueSelectors.\n";
		oss << "\"" << input_string << "\"\n";
		oss << "Error message from Parser:\n";
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	ResidueSelectorOP selector;
	try {
		selector = create_selector( expr_ast );
	} catch ( utility::excn::Exception & e ) {
		std::ostringstream oss;
		oss << "Failed to create a residue selector from the input string:\n";
		oss << "\"" << input_string << "\"\n";
		oss << "Error message from the SelectorLogicParser:\n";
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	return selector;
}


}
}
}

