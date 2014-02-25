// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/pose/selection.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <ObjexxFCL/string.functions.hh>
#include <iostream>


void*
my_main( void* )
{
  utility::vector1< int >  res_vector = utility::tools::make_vector1( -5, -4, -3, 1, 2, 3, 1, 2);
  utility::vector1< char > chain_vector = utility::tools::make_vector1( 'A','A','A','A','A','A',' ','B' );
  std::string tag = make_tag_with_dashes( res_vector, chain_vector );
  std::cout << "CHECK IT: " << tag << std::endl;

  bool string_is_ok;
  std::cout << utility::get_resnum_and_chain( tag, string_is_ok ).first << std::endl;
  std::cout << utility::get_resnum_and_chain( tag, string_is_ok ).second << std::endl;

  exit( 0 );
}


int
main( int argc, char * argv [] )
{

  try {


    ////////////////////////////////////////////////////////////////////////////
    // setup
    ////////////////////////////////////////////////////////////////////////////
    devel::init(argc, argv);

    ////////////////////////////////////////////////////////////////////////////
    // end of setup
    ////////////////////////////////////////////////////////////////////////////

    protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
  }

}
