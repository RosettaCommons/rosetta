// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/Graph.cxxtest.hh
/// @brief  test suite for core::graph::Graph.cc
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/expression_parser/Arithmetic.hh>

/// C++ headers

//Auto Headers


using namespace numeric::expression_parser;


// --------------- Test Class --------------- //

class ArithmeticTests : public CxxTest::TestSuite {

	public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// --------------- Test Cases --------------- //

	void test_parse_simple() {
		std::string test_string1( "1 + 3");
		std::string test_string2( "1 + 3 * 2");
		std::string test_string3( "1 * 3 + 2");
		std::string test_string4( "1 * ( 3 + 2 )");
		ArithmeticScanner as;
		TokenSetOP tokens1 = as.scan( test_string1 );
		TokenSetOP tokens2 = as.scan( test_string2 );
		TokenSetOP tokens3 = as.scan( test_string3 );
		TokenSetOP tokens4 = as.scan( test_string4 );
		std::string tokstr1 = tokens1->print();
		std::string tokstr2 = tokens2->print();
		std::string tokstr3 = tokens3->print();
		std::string tokstr4 = tokens4->print();

		std::string gold_tokstr1 = "Tokens:\nLITERAL(1)\nPLUS_SYMBOL\nLITERAL(3)\n";
		std::string gold_tokstr2 = "Tokens:\nLITERAL(1)\nPLUS_SYMBOL\nLITERAL(3)\nMULTIPLY_SYMBOL\nLITERAL(2)\n";
		std::string gold_tokstr3 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nLITERAL(3)\nPLUS_SYMBOL\nLITERAL(2)\n";
		std::string gold_tokstr4 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nLEFT_PAREN\nLITERAL(3)\nPLUS_SYMBOL\nLITERAL(2)\nRIGHT_PAREN\n";

		TS_ASSERT( gold_tokstr1 == tokstr1 );
		TS_ASSERT( gold_tokstr2 == tokstr2 );
		TS_ASSERT( gold_tokstr3 == tokstr3 );
		TS_ASSERT( gold_tokstr4 == tokstr4 );

		//std::cout << "1:\n" << tokstr1 << std::endl;
		//std::cout << "2:\n" << tokstr2 << std::endl;
		//std::cout << "3:\n" << tokstr3 << std::endl;
		//std::cout << "4:\n" << tokstr4 << std::endl;
		ArithmeticASTExpressionOP expr1( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr2( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr3( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr4( new ArithmeticASTExpression );
		expr1->parse( *tokens1 );
		expr2->parse( *tokens2 );
		expr3->parse( *tokens3 );
		expr4->parse( *tokens4 );
		ASTPrinter printer;
		printer.pretty( false );
		std::string ast_string1 = printer.ast_string( *expr1 );
		std::string ast_string2 = printer.ast_string( *expr2 );
		std::string ast_string3 = printer.ast_string( *expr3 );
		std::string ast_string4 = printer.ast_string( *expr4 );
		//std::cout << ast_string1 << std::endl << std::endl;
		//std::cout << ast_string2 << std::endl << std::endl;
		//std::cout << ast_string3 << std::endl << std::endl;
		//std::cout << ast_string4 << std::endl << std::endl;

		std::string gold_ast_string1 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string2 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string3 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string4 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() )";

		TS_ASSERT( gold_ast_string1 == ast_string1 );
		TS_ASSERT( gold_ast_string2 == ast_string2 );
		TS_ASSERT( gold_ast_string3 == ast_string3 );
		TS_ASSERT( gold_ast_string4 == ast_string4 );

	}

	void test_parse_variables()
	{
		std::string variable_name1( "fa_atr" );
		std::string variable_name2( "fa_rep" );
		std::string variable_name3( "fa_sol" );

		std::string test_string1( "fa_atr + 3");
		std::string test_string2( "1 + fa_sol * 2");
		std::string test_string3( "1 * fa_rep + fa_sol");
		std::string test_string4( "1 * ( fa_atr - fa_rep )");

		ArithmeticScanner as;
		as.add_variable( variable_name1 );
		as.add_variable( variable_name2 );
		as.add_variable( variable_name3 );
		TokenSetOP tokens1 = as.scan( test_string1 );
		TokenSetOP tokens2 = as.scan( test_string2 );
		TokenSetOP tokens3 = as.scan( test_string3 );
		TokenSetOP tokens4 = as.scan( test_string4 );
		std::string tokstr1 = tokens1->print();
		std::string tokstr2 = tokens2->print();
		std::string tokstr3 = tokens3->print();
		std::string tokstr4 = tokens4->print();

		std::string gold_tokstr1 = "Tokens:\nVARIABLE(fa_atr)\nPLUS_SYMBOL\nLITERAL(3)\n";
		std::string gold_tokstr2 = "Tokens:\nLITERAL(1)\nPLUS_SYMBOL\nVARIABLE(fa_sol)\nMULTIPLY_SYMBOL\nLITERAL(2)\n";
		std::string gold_tokstr3 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nVARIABLE(fa_rep)\nPLUS_SYMBOL\nVARIABLE(fa_sol)\n";
		std::string gold_tokstr4 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nLEFT_PAREN\nVARIABLE(fa_atr)\nSUBTRACT_SYMBOL\nVARIABLE(fa_rep)\nRIGHT_PAREN\n";

		TS_ASSERT( gold_tokstr1 == tokstr1 );
		TS_ASSERT( gold_tokstr2 == tokstr2 );
		TS_ASSERT( gold_tokstr3 == tokstr3 );
		TS_ASSERT( gold_tokstr4 == tokstr4 );

		//std::cout << "1:\n" << tokstr1 << std::endl;
		//std::cout << "2:\n" << tokstr2 << std::endl;
		//std::cout << "3:\n" << tokstr3 << std::endl;
		//std::cout << "4:\n" << tokstr4 << std::endl;

		ArithmeticASTExpressionOP expr1( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr2( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr3( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr4( new ArithmeticASTExpression );
		expr1->parse( *tokens1 );
		expr2->parse( *tokens2 );
		expr3->parse( *tokens3 );
		expr4->parse( *tokens4 );
		ASTPrinter printer;
		printer.pretty( false );
		std::string ast_string1 = printer.ast_string( *expr1 );
		std::string ast_string2 = printer.ast_string( *expr2 );
		std::string ast_string3 = printer.ast_string( *expr3 );
		std::string ast_string4 = printer.ast_string( *expr4 );

		//std::cout << ast_string1 << std::endl << std::endl;
		//std::cout << ast_string2 << std::endl << std::endl;
		//std::cout << ast_string3 << std::endl << std::endl;
		//std::cout << ast_string4 << std::endl << std::endl;

		std::string gold_ast_string1 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_atr )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string2 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_sol )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string3 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_rep )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_sol )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) )";
		std::string gold_ast_string4 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_atr )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:SUBTRACT_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:variable:fa_rep )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() )";

		TS_ASSERT( gold_ast_string1 == ast_string1 );
		TS_ASSERT( gold_ast_string2 == ast_string2 );
		TS_ASSERT( gold_ast_string3 == ast_string3 );
		TS_ASSERT( gold_ast_string4 == ast_string4 );

	}

	void test_parse_functions()
	{
		std::string test_string1( "max( 1, 2 + 3 )");
		std::string test_string2( "max( 0.0, min( 1.0, 3*4 ))");
		std::string test_string3( "1 * sqrt( 3.0 )");
		std::string test_string4( "1 * ( sqrt( 2.0 ) + 4.2 )");

		ArithmeticScanner as;
		TokenSetOP tokens1 = as.scan( test_string1 );
		TokenSetOP tokens2 = as.scan( test_string2 );
		TokenSetOP tokens3 = as.scan( test_string3 );
		TokenSetOP tokens4 = as.scan( test_string4 );
		std::string tokstr1 = tokens1->print();
		std::string tokstr2 = tokens2->print();
		std::string tokstr3 = tokens3->print();
		std::string tokstr4 = tokens4->print();

		std::string gold_tokstr1 = "Tokens:\nFUNCTION(max,2)\nLEFT_PAREN\nLITERAL(1)\nCOMMA\nLITERAL(2)\nPLUS_SYMBOL\nLITERAL(3)\nRIGHT_PAREN\n";
		std::string gold_tokstr2 = "Tokens:\nFUNCTION(max,2)\nLEFT_PAREN\nLITERAL(0)\nCOMMA\nFUNCTION(min,2)\nLEFT_PAREN\nLITERAL(1)\nCOMMA\nLITERAL(3)\nMULTIPLY_SYMBOL\nLITERAL(4)\nRIGHT_PAREN\nRIGHT_PAREN\n";
		std::string gold_tokstr3 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nFUNCTION(sqrt,1)\nLEFT_PAREN\nLITERAL(3)\nRIGHT_PAREN\n";
		std::string gold_tokstr4 = "Tokens:\nLITERAL(1)\nMULTIPLY_SYMBOL\nLEFT_PAREN\nFUNCTION(sqrt,1)\nLEFT_PAREN\nLITERAL(2)\nRIGHT_PAREN\nPLUS_SYMBOL\nLITERAL(4.2)\nRIGHT_PAREN\n";

		TS_ASSERT( gold_tokstr1 == tokstr1 );
		TS_ASSERT( gold_tokstr2 == tokstr2 );
		TS_ASSERT( gold_tokstr3 == tokstr3 );
		TS_ASSERT( gold_tokstr4 == tokstr4 );

		//std::cout << "1:\n" << tokstr1 << std::endl;
		//std::cout << "2:\n" << tokstr2 << std::endl;
		//std::cout << "3:\n" << tokstr3 << std::endl;
		//std::cout << "4:\n" << tokstr4 << std::endl;

		ArithmeticASTExpressionOP expr1( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr2( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr3( new ArithmeticASTExpression );
		ArithmeticASTExpressionOP expr4( new ArithmeticASTExpression );
		expr1->parse( *tokens1 );
		expr2->parse( *tokens2 );
		expr3->parse( *tokens3 );
		expr4->parse( *tokens4 );
		ASTPrinter printer;
		printer.pretty( false );
		std::string ast_string1 = printer.ast_string( *expr1 );
		std::string ast_string2 = printer.ast_string( *expr2 );
		std::string ast_string3 = printer.ast_string( *expr3 );
		std::string ast_string4 = printer.ast_string( *expr4 );

		///std::cout << ast_string1 << std::endl << std::endl;
		///std::cout << ast_string2 << std::endl << std::endl;
		///std::cout << ast_string3 << std::endl << std::endl;
		//std::cout << ast_string4 << std::endl << std::endl;

		std::string gold_ast_string1 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTFunction:max:2( ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ), ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) ) )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() )";
		std::string gold_ast_string2 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTFunction:max:2( ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:0 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ), ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTFunction:min:2( ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ), ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTValue:literal:4 )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() )";
		std::string gold_ast_string3 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTFunction:sqrt:1( ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:3 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() )";
		std::string gold_ast_string4 = "ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:1 )  ArithmeticASTRestTerm:MULTIPLY_SYMBOL( ArithmeticASTFactor(  ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTFunction:sqrt:1( ArithmeticASTExpression( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:2 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression:PLUS_SYMBOL( ArithmeticASTTerm( ArithmeticASTFactor(  ArithmeticASTValue:literal:4.2 )  ArithmeticASTRestTerm() ) ArithmeticASTRestExpression() ) ) )  ArithmeticASTRestTerm() ) ) ArithmeticASTRestExpression() )";

		TS_ASSERT( gold_ast_string1 == ast_string1 );
		TS_ASSERT( gold_ast_string2 == ast_string2 );
		TS_ASSERT( gold_ast_string3 == ast_string3 );
		TS_ASSERT( gold_ast_string4 == ast_string4 );

	}

	void test_my_addition() {
		ExpressionCOP add12 = parse_string_to_expression( "1+2" );
		TS_ASSERT( 3 == (*add12)() );
	}

	void test_my_subtraction() {
		ExpressionCOP subtract43 = parse_string_to_expression( "4 - 3" );
		TS_ASSERT( 1 == (*subtract43)() );
	}

	void test_my_multiplication() {
		ExpressionCOP multiply33 = parse_string_to_expression( "3 * 3" );
		TS_ASSERT( 9 == (*multiply33)() );
	}

	void test_my_division() {
		ExpressionCOP divide24 = parse_string_to_expression( "2 / 4" );
		TS_ASSERT( 0.5 == (*divide24)() );

	}

	void test_multiple_additions() {
		ExpressionCOP add123 = parse_string_to_expression( "1+2+3" );
		TS_ASSERT( 6 == (*add123)() );
	}

	void test_addition_order() {
		/// If this were parsed as 12 - (5+3) then the answer would be 4 instead of 10.
		ExpressionCOP sub1253 = parse_string_to_expression( "12 - 5 + 3" );
		TS_ASSERT( 10 == (*sub1253)() );
	}

	void test_multiplication_order() {
		/// If this were parsed as 16 / (4*2) then the answer would be 2 instead of 8.
		ExpressionCOP divmult = parse_string_to_expression( "16 / 4 * 2" );
		TS_ASSERT( 8 == (*divmult)() );
	}

	void test_precedence() {
 		ExpressionCOP precidence = parse_string_to_expression( "12 + 5 * 3" );
		TS_ASSERT( 27 == (*precidence)() );
	}

	void test_parens() {
 		ExpressionCOP parens = parse_string_to_expression( "(12 + 5) * 3" );
		TS_ASSERT( 51 == (*parens)() );
	}

	void test_min() {
		ExpressionCOP minex = parse_string_to_expression( "min( 12, 5 )" );
		TS_ASSERT( 5 == (*minex)() );
	}

	void test_max() {
		ExpressionCOP maxex = parse_string_to_expression( "max( 12, 5 )" );
		TS_ASSERT( 12 == (*maxex)() );
	}

	void test_sqrt() {
		ExpressionCOP sqrtex = parse_string_to_expression( "sqrt( 36 )" );
		TS_ASSERT( 6 == (*sqrtex)() );
		sqrtex = parse_string_to_expression( "sqrt( 4 )" );
		TS_ASSERT( 2 == (*sqrtex)() );
		sqrtex = parse_string_to_expression( "sqrt( 25 )" );
		TS_ASSERT( 5 == (*sqrtex)() );
	}

	void test_combine_functions() {
		ExpressionCOP boundex = parse_string_to_expression( "max( 0.0, min( 1.0, 5 ) )" );
		TS_ASSERT( 1.0 == (*boundex)() );
		boundex = parse_string_to_expression( "max( 0.0, min( 1.0, 0.5 ) )" );
		TS_ASSERT( 0.5 == (*boundex)() );
		boundex = parse_string_to_expression( "max( 0.0, min( 1.0, -12 ) )" );
		TS_ASSERT( 0.0 == (*boundex)() );
	}

	void test_mix_addition_and_functions() {
		ExpressionCOP combo = parse_string_to_expression( "max( 1, 3 ) + 3"  );
		TS_ASSERT( 6 == (*combo)() );
		combo = parse_string_to_expression( "max( 1+3, 2 ) + 1" );
		TS_ASSERT( 5 == (*combo)() );
	}

	void test_mix_multiplication_and_functions()
	{
		ExpressionCOP combo = parse_string_to_expression( "max( 1, 3 ) * 4" );
		TS_ASSERT( 12 == (*combo)() );
		combo = parse_string_to_expression( "4 * max( 1, 3 )" );
		TS_ASSERT( 12 == (*combo)() );
		combo = parse_string_to_expression( "4 * max( 2*2, 3*3 )" );
		TS_ASSERT( 36 == (*combo)() );
	}

	void test_function_order() {
		ExpressionCOP order = parse_string_to_expression( "1 - max( 1, 3 ) + 5" );
		TS_ASSERT( 3 == (*order)() );
		order = parse_string_to_expression( "5 / max( 2, 4 ) * 8" );
		TS_ASSERT( 10 == (*order)() );
	}

	void test_function_precidence() {
		ExpressionCOP precidence = parse_string_to_expression( "1 + max(1,3)*5" );
		TS_ASSERT( 16 == (*precidence)() );
		precidence = parse_string_to_expression( "max(1,3)*5 + 1 " );
		TS_ASSERT( 16 == (*precidence)() );
		precidence = parse_string_to_expression( "min(14,16) + 5*3" );
		TS_ASSERT( 29 == (*precidence)() );
	}


	ExpressionCOP
	parse_string_with_variables_to_expression(
		std::string const & input_string,
		SimpleExpressionCreator & sec
	) const
	{
		ArithmeticScanner as;
		std::map< std::string, VariableExpressionOP > namemap = sec.variables();
		for ( std::map< std::string, VariableExpressionOP >::const_iterator
				iter = namemap.begin(), iter_end = namemap.end();
				iter != iter_end; ++iter ) {
			as.add_variable( iter->first );
		}
		TokenSetOP tokens = as.scan( input_string );
		ArithmeticASTExpression expr_ast;
		expr_ast.parse( *tokens );
		return sec.create_expression_tree( expr_ast );
	}


	std::list< std::string >
	regular_varnames() {
		std::list< std::string > return_list;
		return_list.push_back( "x" );
		return_list.push_back( "y" );
		return_list.push_back( "z" );
		return return_list;
	}

	void test_add_derivative() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP addxy = parse_string_with_variables_to_expression( "x + y", sec );
		sec.get_variable( "x" )->set_value( 2 );
		sec.get_variable( "y" )->set_value( 3 );
		ExpressionCOP daddxy_dx = addxy->differentiate( "x" );
		TS_ASSERT( daddxy_dx );
		TS_ASSERT( 1 == (*daddxy_dx)() );
		ExpressionCOP daddxy_dy = addxy->differentiate( "y" );
		TS_ASSERT( daddxy_dy );
		TS_ASSERT( 1 == (*daddxy_dy)() );
	}

	void test_multiply_derivative() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP mul3x = parse_string_with_variables_to_expression( "3 * x", sec );
		sec.get_variable( "x" )->set_value( 2 );
		ExpressionCOP mul3x_dx = mul3x->differentiate( "x" );
		TS_ASSERT( 3 == (*mul3x_dx)() );
		ExpressionCOP mul3x_dy = mul3x->differentiate( "y" );
		TS_ASSERT( mul3x_dy == 0 );
	}

	void test_multiply_and_add_derivatives() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP addmul = parse_string_with_variables_to_expression( "3 * x + 5 * y", sec );
		sec.get_variable( "x" )->set_value( 2 );
		sec.get_variable( "y" )->set_value( 7 );
		ExpressionCOP addmul_dx = addmul->differentiate( "x" );
		TS_ASSERT( 3 == (*addmul_dx)() );
		ExpressionCOP addmul_dy = addmul->differentiate( "y" );
		TS_ASSERT( 5 == (*addmul_dy)() );
	}

	void test_multiply_and_add_derivatives2() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP addmul = parse_string_with_variables_to_expression( "3 * x * y + 5 * y * z", sec );
		sec.get_variable( "x" )->set_value( 2 );
		sec.get_variable( "y" )->set_value( 7 );
		sec.get_variable( "z" )->set_value( 11 );
		ExpressionCOP addmul_dx = addmul->differentiate( "x" );
		TS_ASSERT( 21 == (*addmul_dx)() );
		ExpressionCOP addmul_dy = addmul->differentiate( "y" );
		TS_ASSERT( 61 == (*addmul_dy)() );
		ExpressionCOP addmul_dz = addmul->differentiate( "z" );
		TS_ASSERT( 35 == (*addmul_dz)() );
	}

	void test_equals_expression() {
		ExpressionCOP fourexp = parse_string_to_boolean_expression( "EQUALS( 4, 2 * 2 )" );
		TS_ASSERT( (*fourexp)() == 1.0 );
		ExpressionCOP fiveexp = parse_string_to_boolean_expression( "EQUALS( 4, 2 * 2 + 1 )" );
		TS_ASSERT( (*fiveexp)() == 0.0 );
	}

	void test_gt_expression() {
		ExpressionCOP fourexp = parse_string_to_boolean_expression( "GT( 3, 2 * 2 )" );
		TS_ASSERT( (*fourexp)() == 0.0 );
		ExpressionCOP fiveexp = parse_string_to_boolean_expression( "GT( 6, 2 * 2 + 1 )" );
		TS_ASSERT( (*fiveexp)() == 1.0 );
		ExpressionCOP sixexp = parse_string_to_boolean_expression( "GT( 6, 2 * ( 2 + 1 ) )" );
		TS_ASSERT( (*sixexp)() == 0.0 );
	}

	void test_gte_expression() {
		ExpressionCOP fourexp = parse_string_to_boolean_expression( "GTE( 3, 2 * 2 )" );
		TS_ASSERT( (*fourexp)() == 0.0 );
		ExpressionCOP fiveexp = parse_string_to_boolean_expression( "GTE( 6, 2 * 2 + 1 )" );
		TS_ASSERT( (*fiveexp)() == 1.0 );
		ExpressionCOP sixexp = parse_string_to_boolean_expression( "GTE( 6, 2 * ( 2 + 1 ) )" );
		TS_ASSERT( (*sixexp)() == 1.0 );
	}

	void test_lt_expression() {
		ExpressionCOP fourexp = parse_string_to_boolean_expression( "LT( 3, 2 * 2 )" );
		TS_ASSERT( (*fourexp)() == 1.0 );
		ExpressionCOP fiveexp = parse_string_to_boolean_expression( "LT( 6, 2 * 2 + 1 )" );
		TS_ASSERT( (*fiveexp)() == 0.0 );
		ExpressionCOP sixexp = parse_string_to_boolean_expression( "LT( 6, 2 * ( 2 + 1 ) )" );
		TS_ASSERT( (*sixexp)() == 0.0 );
	}

	void test_lte_expression() {
		ExpressionCOP fourexp = parse_string_to_boolean_expression( "LTE( 3, 2 * 2 )" );
		TS_ASSERT( (*fourexp)() == 1.0 );
		ExpressionCOP fiveexp = parse_string_to_boolean_expression( "LTE( 6, 2 * 2 + 1 )" );
		TS_ASSERT( (*fiveexp)() == 0.0 );
		ExpressionCOP sixexp = parse_string_to_boolean_expression( "LTE( 6, 2 * ( 2 + 1 ) )" );
		TS_ASSERT( (*sixexp)() == 1.0 );
	}

	void test_and_expession() {
		ExpressionCOP zerozero = parse_string_to_boolean_expression( "AND( EQUALS( 1.0, 0.5 ), EQUALS( 4, 2 * 2 + 1 ))" );
		ExpressionCOP zeroone  = parse_string_to_boolean_expression( "AND( EQUALS( 1.0, 0.5 ), EQUALS( 4, 2 * 2 ))" );
		ExpressionCOP onezero  = parse_string_to_boolean_expression( "AND( EQUALS( 1.0, 1.0 ), EQUALS( 4, 2 * 2 + 1 ))" );
		ExpressionCOP oneone   = parse_string_to_boolean_expression( "AND( EQUALS( 1.0, 1.0 ), EQUALS( 4, 2 * 2 ))" );

		TS_ASSERT( (*zerozero)() == 0.0 );
		TS_ASSERT( (*zeroone)()  == 0.0 );
		TS_ASSERT( (*onezero)()  == 0.0 );
		TS_ASSERT( (*oneone)()   == 1.0 );
	}

	void test_or_expession() {
		ExpressionCOP zerozero = parse_string_to_boolean_expression( "OR( EQUALS( 1.0, 0.5 ), EQUALS( 4, 2 * 2 + 1 ))" );
		ExpressionCOP zeroone  = parse_string_to_boolean_expression( "OR( EQUALS( 1.0, 0.5 ), EQUALS( 4, 2 * 2 ))" );
		ExpressionCOP onezero  = parse_string_to_boolean_expression( "OR( EQUALS( 1.0, 1.0 ), EQUALS( 4, 2 * 2 + 1 ))" );
		ExpressionCOP oneone   = parse_string_to_boolean_expression( "OR( EQUALS( 1.0, 1.0 ), EQUALS( 4, 2 * 2 ))" );

		TS_ASSERT( (*zerozero)() == 0.0 );
		TS_ASSERT( (*zeroone)()  == 1.0 );
		TS_ASSERT( (*onezero)()  == 1.0 );
		TS_ASSERT( (*oneone)()   == 1.0 );
	}

	void test_active_variables_simple() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP xplusy = parse_string_with_variables_to_expression( "x + y", sec );

		std::list< std::string > active = xplusy->active_variables();
		TS_ASSERT( find(active.begin(), active.end(), "x") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") != active.end() )


		ExpressionCOP sqrtx = parse_string_with_variables_to_expression( "sqrt(x)", sec );

		active = sqrtx->active_variables();
		TS_ASSERT( find(active.begin(), active.end(), "x") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") == active.end() )


	}

	void test_active_variables_max() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP max_xy = parse_string_with_variables_to_expression( "max( x, y )", sec );
		sec.get_variable("x")->set_value( 3 );
		sec.get_variable("y")->set_value( 2 );
		std::list< std::string > active = max_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") == active.end() )

		sec.get_variable("x")->set_value( 30 );
		sec.get_variable("y")->set_value( 32 );
		active = max_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") == active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") != active.end() )

	}

	void test_active_variables_min() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP min_xy = parse_string_with_variables_to_expression( "min( x, y )", sec );
		sec.get_variable("x")->set_value( 3 );
		sec.get_variable("y")->set_value( 2 );
		std::list< std::string > active = min_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") == active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") != active.end() )

		sec.get_variable("x")->set_value( 30 );
		sec.get_variable("y")->set_value( 32 );
		active = min_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") == active.end() )

	}


	void test_active_variables_min_complicated() {
		std::list< std::string > vars( regular_varnames() );
		SimpleExpressionCreator sec( vars );
		ExpressionCOP min_xy = parse_string_with_variables_to_expression( "min( x+z, y )", sec );
		sec.get_variable("x")->set_value( 3 );
		sec.get_variable("y")->set_value( 2 );
		sec.get_variable("z")->set_value( -2 );
		std::list< std::string > active = min_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") == active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "z") != active.end() )

		sec.get_variable("x")->set_value( 30 );
		sec.get_variable("y")->set_value( 32 );
		sec.get_variable("z")->set_value( 4 );
		active = min_xy->active_variables();

		TS_ASSERT( find(active.begin(), active.end(), "x") == active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "y") != active.end() )
		TS_ASSERT( find(active.begin(), active.end(), "z") == active.end() )

	}


};


