// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// Based on the Boost::Spirit calculator examples.
// (See also rosetta/rosetta_source/external/boost_*/LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)

/// @file   core/fragment/picking/Calculator.cc
/// @brief a string-input calculator
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#include <numeric/Calculator.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/util.hh>

#include <utility/vector1.hh>
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>

namespace numeric {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::qi::ascii;
namespace phoenix = boost::phoenix;

class CalculatorParser;

// This strange, indirect invocation is to make boost::phoenix::bind happy.
// If you can figure out how to directly call the member function, have at it.
void do_add_symbol(CalculatorParser & cp, std::string name, double value);

// These adapter functions are needed to get the type determination for the underlying overloaded functions correct
double do_abs( double a) { return std::abs(a); }
double do_pow( double a, double b ) { return std::pow(a, b); }
double do_exp( double a ) { return std::exp(a); }
double do_ln( double a ) { return std::log(a); }
double do_log10( double a ) { return std::log10(a); }
double do_log2( double a ) { return numeric::log(a,2.0); }
double do_log( double a, double b ) { return  numeric::log(a,b); }
double do_sqrt( double a ) { return std::sqrt(a); }
double do_sin( double a ) { return std::sin(a); }
double do_cos( double a ) { return std::cos(a); }
double do_tan( double a ) { return std::tan(a); }
double do_max( std::vector<double> a ) { return numeric::max( utility::vector1<double>(a) ); }
double do_min( std::vector<double> a ) { return numeric::min( utility::vector1<double>(a) ); }
double do_mean( std::vector<double> a ) { return numeric::mean( utility::vector1<double>(a) ); }
double do_median( std::vector<double> a) { return numeric::median( utility::vector1<double>(a) ); }

class CalculatorParser : qi::grammar<std::string::iterator, double(), ascii::space_type> {

public:
	typedef qi::rule<std::string::iterator, double(), ascii::space_type> Rule;
	typedef qi::rule<std::string::iterator, std::string(), ascii::space_type> RuleS;
	typedef qi::rule<std::string::iterator, void(), ascii::space_type> RuleX;

	typedef std::map<std::string,numeric::Real> Varmap;

	// Setup the grammar of the calculator
	CalculatorParser(Varmap variables) :
		CalculatorParser::base_type(parser_)
	{

		for ( Varmap::iterator iter(variables.begin()); iter != variables.end(); ++iter ) {
			add_symbol(iter->first,iter->second);
		}

		using namespace boost::spirit::qi;
		using boost::spirit::qi::_1;

		varname = as_string[ lexeme[ alpha >> *( alnum | '_' ) ]];

		assignment = omit[ ( varname >> '=' >> expression ) [ phoenix::bind(do_add_symbol, boost::ref(*this), qi::_1, qi::_2) ] ];

		expression =
			term                            [_val = _1]
			>> *(   ('+' >> term            [_val += _1])
			|   ('-' >> term            [_val -= _1])
		)
			;

		term =
			exponent                        [_val = _1]
			>> *(   ('*' >> exponent        [_val *= _1])
			|   ('/' >> exponent        [_val /= _1])
		)
			;

		exponent =
			factor                        [_val = _1]
			>> -('^' >> exponent   [ _val = phoenix::bind(do_pow, qi::_val, qi::_1) ] )
			;

		factor =
			double_                         [_val = _1]
			| function                      [_val = _1]
			| ( local_vars_ >> ! lit('(') ) [_val = _1]
			|   '(' >> expression           [_val = _1] >> ')'
			|   ('-' >> factor              [_val = -_1])
			|   ('+' >> factor              [_val = _1])
			;

		function =
			( no_case["abs"] >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_abs, qi::_1) ]
			| ( no_case["exp"] >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_exp, qi::_1) ]
			| ( no_case["ln"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_ln, qi::_1) ]
			| ( no_case["log"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_log10, qi::_1) ]
			| ( no_case["log"]  >> '(' >> expression >> ',' >> expression >> ')' ) [ _val = phoenix::bind(do_log, qi::_1, qi::_2) ]
			| ( no_case["log10"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_log10, qi::_1) ]
			| ( no_case["log2"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_log2, qi::_1) ]
			| ( no_case["r2d"]  >> '(' >> expression >> ')' ) [ _val = _1 * numeric::NumericTraits< double >::rad2deg() ]
			| ( no_case["d2r"]  >> '(' >> expression >> ')' ) [ _val = _1 * numeric::NumericTraits< double >::deg2rad() ]
			| ( no_case["sqrt"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_sqrt, qi::_1) ]
			| ( no_case["cos"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_cos, qi::_1) ]
			| ( no_case["sin"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_sin, qi::_1) ]
			| ( no_case["tan"]  >> '(' >> expression >> ')' ) [ _val = phoenix::bind(do_tan, qi::_1) ]
			// indeterminite argument functions
			| ( no_case["mean"]  >> '(' >> (expression % ',') >> ')' ) [ _val = phoenix::bind(do_mean, qi::_1) ]
			| ( no_case["median"]  >> '(' >> (expression % ',') >> ')' ) [ _val = phoenix::bind(do_median, qi::_1) ]
			| ( no_case["min"]  >> '(' >> (expression % ',') >> ')' ) [ _val = phoenix::bind(do_min, qi::_1) ]
			| ( no_case["max"]  >> '(' >> (expression % ',') >> ')' ) [ _val = phoenix::bind(do_max, qi::_1) ]
			;

		parser_ = *(assignment >> ';') >> expression; //assigment is omitted, so expression is the value for all.
	}

	void add_symbol(std::string name, double value) {
		local_vars_.add(name, value);
	}

	Rule parser() { return parser_; }

private:
	RuleS varname;
	RuleX assignment;
	Rule expression, term, exponent, factor, function;
	Rule parser_;
	qi::symbols<char, double> local_vars_;
};

void do_add_symbol(CalculatorParser & cp, std::string name, double value) {
	cp.add_symbol(name,value);
}

//--------- Functions for externally visible calculator object ----------------///

Calculator::~Calculator() {}

Calculator::Calculator(std::string const & equation):
	equation_(equation)
{}

// Return true on failure and false on success
bool
Calculator::compute( std::map<std::string, Real> & values, Real & output) const {
	std::string equation( equation_ );
	std::string::iterator iter = equation.begin();
	CalculatorParser calc_parser( values );

	bool success = qi::phrase_parse(iter, equation.end(), calc_parser.parser(), ascii::space, output);
	if ( success && iter == equation.end() ) {
		return false;
	}
	//std::cout << "Output -- " << output << " and parsing " << (success?"Succeeded":"Failed") << std::endl;
	//std::cout << "'" << std::string(equation.begin(), iter) << "' and '" << std::string(iter, equation.end()) << "'" << std::endl;
	return true;

}

} // numeric

