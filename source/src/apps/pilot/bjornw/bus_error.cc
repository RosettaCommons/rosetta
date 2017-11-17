// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include <utility/io/ozstream.hh>
#include <iostream>
#include <fstream>

int
main( int argc, char * argv [] )
{
  try{
    utility::io::ozstream(output);
    output.open("bus_error3");
    output << "TEST1\n"  << std::endl;
    output.close();

    std::ofstream output2 ("bus_error.std3");
    output2 << "TEST2\n";
    output2.close();

  } catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }

  return 0;
}
