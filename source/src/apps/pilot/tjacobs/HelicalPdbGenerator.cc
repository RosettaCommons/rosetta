// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelicalPdbGenerator.cc
///
/// @brief

/// @author Tim jacobs

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//Core
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

//C++ headers
#include <iostream>
#include <string>
#include <stdlib.h>

class HelicalPdbGeneratorMover : public protocols::moves::Mover {

public:

  HelicalPdbGeneratorMover();

  virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
    return new HelicalPdbGeneratorMover( *this );
  }

  virtual
  std::string
  get_name() const {
    return "HelicalPdbGeneratorMover";
  }

  virtual protocols::moves::MoverOP fresh_instance() const {
    return new HelicalPdbGeneratorMover();
  }

};

HelicalPdbGeneratorMover::HelicalPdbGeneratorMover() {}

void HelicalPdbGeneratorMover::apply(core::pose::Pose & pose){
  //set up dssp info - necessary in order to find helices based on secondary structure
	//....duh

  if(pose.total_residue() > 0){

    core::scoring::dssp::Dssp dssp( pose );
    dssp.insert_ss_into_pose( pose );

    utility::vector1< std::pair< Size,Size > > helix_endpts;
    for(Size i=1; i<=pose.total_residue(); i++){

          //find all the strands in the structure
          Size helix_start(0);
          Size helix_end(0);
          if(pose.secstruct(i) == 'H'){

              helix_start=i;
              while(i<=pose.total_residue() && pose.secstruct(i)=='H'){
                  i++;
              }
              helix_end=i;

              helix_endpts.push_back(std::make_pair(helix_start, helix_end));
          }
      }

    Pose outputPose;
    for(Size i=1; i<=helix_endpts.size(); ++i){
        outputPose.append_residue_by_jump(pose.residue(helix_endpts[i].first), outputPose.total_residue(), "", "", false);
        for(Size j=helix_endpts[i].first; j<=helix_endpts[i].second; ++j){
            outputPose.append_residue_by_bond(pose.residue(j));
        }
    }

    pose = outputPose;
  }
}

int
main( int argc, char * argv [] )
{

	try {

  using namespace std;
  using namespace utility;
  using namespace core;

  // initialize core
  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go( new HelicalPdbGeneratorMover() );

  std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
