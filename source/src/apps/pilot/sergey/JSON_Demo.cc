// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief  Demo for JSON parser
/// @author Sergey Lyskov


#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>

#include <utility/tools/make_vector.hh>
#include <utility/tools/make_map.hh>

#include <string>
#include <iostream>

#include <utility/excn/Exceptions.hh>

int main( int argc, char * argv [] )
{

	try {

	utility::json_spirit::Value v;

	std::string input = "5.55";

	std::cout << "Parsing: " << input << std::endl;
	std::cout << "Result:" << utility::json_spirit::read(input, v) << std::endl;
	std::cout << "Result Type:" << v.type() << " get_real="<< v.get_real() << std::endl;
	std::cout << std::endl << std::endl;


	input = " <XML> ";
	std::cout << "Parsing: " << input << std::endl;
	std::cout << "Result:" << utility::json_spirit::read(input, v) << std::endl;
	std::cout << std::endl << std::endl;


	input = " { \"A\" : [1, 2, 3] } ";
	utility::json_spirit::mValue mv;
	std::cout << "Parsing: " << input << std::endl;
	std::cout << "Result:" << utility::json_spirit::read(input, mv) << std::endl;
	std::cout << "RootType=" << mv.type() << std::endl;
	if(mv.get_obj().count("A")) {
		std::cout << "Key A exist!!!" << std::endl;
		std::cout << "mv.get_obj()[A]: [";

		utility::json_spirit::mArray const & array = mv.get_obj()["A"].get_array();

		for(int i; i<array.size(); ++i) {
			std::cout << array[i].get_int() << ", ";
		}
		std::cout << "]" << std::endl;
	}

	std::cout << "JSON Emmitter demo..." << std::endl;
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;
	using utility::tools::make_vector;
	using utility::tools::make_map;

	std::cout << utility::json_spirit::write(make_vector(Value(1), Value(2), Value(3), Value(5)), utility::json_spirit::pretty_print) << std::endl;

	/* Let create a Loop-like object, something like this:
		"Loop" : {
			"Start" : {"ResNo" : "50", "Insertion" : "A", "Chain" : "A"},
			"Cutpoint" : {"ResNo" : "50", "Insertion" : "A", "Chain" : "A"},
			"Stop" : {"ResNo" : "50", "Insertion" : "A", "Chain" : "A"},
		}
	*/
	Value loop( make_vector( Pair("Loop",
								          make_vector( Pair("Start",    make_vector(Pair("ResNo", 50), Pair("Insertion", "A"), Pair("Chain", "A") ) )
													  ,Pair("Cutpoint", make_vector(Pair("ResNo", 50), Pair("Insertion", "A"), Pair("Chain", "A") ) )
													  ,Pair("Stop",     make_vector(Pair("ResNo", 50), Pair("Insertion", "A"), Pair("Chain", "A") ) )
								          			  ) ) ) );

	std::cout << utility::json_spirit::write(loop, utility::json_spirit::pretty_print) << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

