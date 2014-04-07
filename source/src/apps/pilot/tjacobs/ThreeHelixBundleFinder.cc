//
//  ThreeHelixBundleFinder.cc
//  Rosetta
//
//  Created by Tim Jacobs on 1/23/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <devel/init.hh>
#include <devel/sewing/ThreeHelixBundleFinderMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

int
main( int argc, char * argv [] )
{

  try {

    devel::init(argc, argv);
    protocols::jd2::JobDistributor::get_instance()->go( new ThreeHelixBundleFinderMover() );

  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
    return -1;
  }

}
