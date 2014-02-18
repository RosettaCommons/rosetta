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

#include <protocols/jobdist/standard_mains.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/init/init.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/util.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <protocols/jd2/JobDistributor.hh>
#include <core/init/init.hh>

//#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


//using basic::T;
//using core::util::Error;
//using core::util::Warning;


//using namespace core;
//using namespace protocols;

//using utility::vector1;



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // init option system
  core::init::init(argc, argv);

  // Make a sequence move to hold our pack rotamers mover and relax mover
  protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;

  //create a task factory: this will create a new PackerTask for each input pose
  core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
  main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
  if ( option[ packing::resfile ].user() ) {
    main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
  }

  //create a ScoreFunction from commandline options
  core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();
  
  //create the PackRotamersMover which will do the packing
  protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
  pack_mover->task_factory( main_task_factory );
  pack_mover->score_function( score_fxn );

  // add the pack rotamers mover to the sequnce mover
  seqmov->add_mover( pack_mover );

  // create the relax mover
  protocols::moves::MoverOP relax_mover = protocols::relax::generate_relax_from_cmd();

  // add the relax move to the sequence mover
  seqmov->add_mover( relax_mover );

  // get and instance of the job distruibutor and pass it the sequence mover
  protocols::jd2::JobDistributor::get_instance()->go(seqmov);

  // Done
  return 0;
}
