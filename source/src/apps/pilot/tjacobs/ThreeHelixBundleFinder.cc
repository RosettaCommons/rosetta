// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Tim Jacobs

#include <devel/init.hh>
#include <protocols/sewing/ThreeHelixBundleFinderMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

int
main( int argc, char * argv [] )
{

  try {

    devel::init(argc, argv);
    protocols::jd2::JobDistributor::get_instance()->go( new ThreeHelixBundleFinderMover() );

  } catch (utility::excn::Exception const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
    return -1;
  }

}
