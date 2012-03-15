// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief  Demo for JSON parser
/// @author Sergey Lyskov


#include <utility/json_spirit/json_spirit_reader.h>

#include <string>
#include <iostream>

int main( int argc, char * argv [] )
{
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


	input = " \"A\" : [1, 2, 3] ";
	//utility::json_spirit::Value mv;
	std::cout << "Parsing: " << input << std::endl;
	std::cout << "Result:" << utility::json_spirit::read(input, v) << std::endl;
	/* if(v.get_obj().count("A")) {
		std::cout << "Key A exist!!!" << v.get_obj()["A"] << std::endl;
	} */

}

