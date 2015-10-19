// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FoldTreeOutputter.cc
///
/// @brief

/// @author Tim jacobs

// Core headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

//C++ headers
#include <iostream>
#include <string>
#include <stdlib.h>

class FoldTreeOutputter : public protocols::moves::Mover {

public:

  FoldTreeOutputter();

  virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
		//Stupid and useless comment
    return new FoldTreeOutputter( *this );
  }

  virtual
  std::string
  get_name() const {
    return "FoldTreeOutputter";
  }

  virtual protocols::moves::MoverOP fresh_instance() const {
    return new FoldTreeOutputter();
  }

};

FoldTreeOutputter::FoldTreeOutputter() {}

void FoldTreeOutputter::apply(core::pose::Pose & pose){
  std::cout << pose.fold_tree() << std::endl;
  std::cout << "------DONE------" << std::endl;
}

int
main( int argc, char * argv [] )
{

	try {

	//Another useless comment
  using namespace std;
  using namespace utility;
  using namespace core;

  // initialize core
  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go( new FoldTreeOutputter() );

  std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
