// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FoldTreeOutputter.cc
///
/// @brief

/// @author Tim jacobs

// Core headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/dssp/Dssp.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

#include <utility/excn/Exceptions.hh>

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

	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );
	std::cout << "Working on job: " << job_name << std::endl;

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->score(pose);

	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();

	std::string secstruct = dssp.get_dssp_secstruct();

	core::Size sum_neighbors = 0;
	core::Size cur_seg_size = 0;
	for(core::Size i=1; i<=pose.size(); ++i) {

		if(cur_seg_size == 0 ){
			cur_seg_size = 0;
		}
		else if (dssp.get_dssp_secstruct(i) != dssp.get_dssp_secstruct(i-1)){

		}

		++cur_seg_size;

		core::Size num_neighbors = pose.energies().tenA_neighbor_graph().get_node(i)->num_neighbors_counting_self();
		std::cout << "Num neighbors for res " << i << ": " << num_neighbors << std::endl;
		sum_neighbors += num_neighbors;
	}
	std::cout << "Num neighbors full protein: " << sum_neighbors << std::endl;
	std::cout << "Average neighbors: " << sum_neighbors/pose.size() << std::endl;


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

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
