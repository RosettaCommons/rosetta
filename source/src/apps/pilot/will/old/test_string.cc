// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

int
main (int /*argc*/, char ** /* *argv[] */)
{

	using utility::to_string;
	using utility::from_string;
	using std::string;

	std::cout << "testing to_string" << std::endl;

	std::cout << "'" << to_string<double>(5.0 ) << "'" << std::endl;
	std::cout << "'" << to_string<float>(5.0f) << "'" << std::endl;
	std::cout << "'" << to_string<int>(5   ) << "'" << std::endl;
	std::cout << "'" << to_string<char>('5' ) << "'" << std::endl;
	std::cout << "'" << to_string<string>("5" ) << "'" << std::endl;

	std::cout << "testing from_string" << std::endl;

	// char c = from_string<char>( std::string('c') );
	// std::cout << "'" << c << "'" << std::endl;

	// float f = from_string<float>(std::string("5.84593"));
	// std::cout << "'" << f << "'" << std::endl;
	// int i   = from_string<int>  (std::string("5"));
	// std::cout << "'" << i << "'" << std::endl;

	return 0;

}
