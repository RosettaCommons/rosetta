// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/select/residue_selector/SelectorLogicParser.hh
/// @brief  Construct selectors as logical combinations of other residue selectors using AND, NOT, and OR
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_SelectorLogicParser_HH
#define INCLUDED_core_select_residue_selector_SelectorLogicParser_HH

// Unit headers
#include <core/select/residue_selector/SelectorLogicParser.fwd.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Numeric headers
#include <numeric/expression_parser/Arithmetic.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace select {
namespace residue_selector {

/// @brief Class to traverse the abstract syntax tree produced
/// by the parsing of a properly-formed string in the Arithmetic
/// expression language.  Produces a ResidueSelector tree capable
/// of logically joining together already-defined ResidueSelectors.
class SelectorLogicParser : public numeric::expression_parser::ASTVisitor
{
public:
	SelectorLogicParser();
	~SelectorLogicParser() override;

	ResidueSelectorOP
	parse_string_to_residue_selector(
		basic::datacache::DataMap const & dm,
		std::string const & selector_logic_string
	);

	ResidueSelectorOP
	parse_string_to_residue_selector(
		std::string const & input_string,
		std::map< std::string, ResidueSelectorOP > const & selectors
	);

	void add_named_selector( std::string const & name, ResidueSelectorOP selector );

	void
	visit( numeric::expression_parser::ArithmeticASTExpression const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTAndClause const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTOrClause const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTFunction const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTTerm const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTFactor const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTValue const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTRestTerm const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTRestAndClause const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTRestOrClause const & ) override;

	void
	visit( numeric::expression_parser::ArithmeticASTRestExpression const & ) override;


	void
	visit( numeric::expression_parser::ArithmeticASTNode const & ) override;

	ResidueSelectorOP
	create_selector( numeric::expression_parser::ArithmeticASTExpression const & );

	/// @brief Factory method to be implemented by derived classes
	/// which may wish to handle variable expressions in a specific manner
	virtual
	ResidueSelectorOP
	handle_variable_expression( numeric::expression_parser::ArithmeticASTValue const & );

	/// @brief Factory method to be implemented by derived classes
	/// which may wish to handle function expressions in a specific manner
	virtual
	ResidueSelectorOP
	handle_function_expression(
		numeric::expression_parser::FunctionTokenCOP function,
		utility::vector1< ResidueSelectorOP > const & args
	);

private:

	/// @brief After creating an ArithmeticScanner and loading both it and
	/// this with the set of available ResidueSelectors, complete the
	/// scanning, parsing, and residue-selector construction steps of converting
	/// a string to a ResidueSelector.
	ResidueSelectorOP
	finish_parse_string_to_residue_selector(
		std::string const & input_string,
		numeric::expression_parser::ArithmeticScanner & as
	);

private:

	// Data

	// The already-delcared residue selectors that will be combined according to the instructions
	// in the AST
	std::map< std::string, ResidueSelectorOP > variables_;

	// Created inside a traversal of the AST -- return sub-selector trees
	// using this variable, and keep living trees alive
	// by storing them on the stack during the AST recursive traversal.
	ResidueSelectorOP last_constructed_selector_;

	// For partially-constructed sub-expressions, use this variable
	ResidueSelectorOP semi_constructed_selector_;

};


}
}
}

#endif
