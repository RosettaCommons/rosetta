// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/stranges/sym_fixbb.cc
/// @brief  Symmetric fixed backbone design

//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <basic/options/option.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
//#include <core/scoring/Energies.hh>


//protocols library (Movers)
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

//symmetry
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <devel/init.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.stranges.sym_fixbb" );

//local options
namespace basic{ namespace options{ namespace OptionKeys{
basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
}}}//basic::options::OptionKeys

using namespace protocols::moves;
class SymFixbbMover : public Mover {
public:
  SymFixbbMover();

  virtual void apply( core::pose::Pose& pose );

  virtual void setup_task_mm(core::pose::Pose & pose, core::pack::task::TaskFactoryOP tf  );

  virtual MoverOP clone() const {
    return new SymFixbbMover( *this );
  }

  virtual MoverOP fresh_instance() const {
    return new SymFixbbMover;
  }

private:
  core::scoring::ScoreFunctionOP scorefxn_a;
	core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_;
  core::kinematics::MoveMapOP movemap_;
	core::pack::task::TaskFactoryOP main_task_factory_;
  core::pack::task::PackerTaskOP task_;
};

//constructor
SymFixbbMover::SymFixbbMover() {
	using namespace basic::options;
  using namespace basic::options::OptionKeys;
	scorefxn_a = core::scoring::get_score_function();
	scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( scorefxn_a );
  movemap_ = new core::kinematics::MoveMap();
	main_task_factory_ = new core::pack::task::TaskFactory;
	main_task_factory_->push_back( new core::pack::task::operation::InitializeFromCommandline );
  if ( option[ packing::resfile ].user() ) {
		TR << "Reading a resfile for packing." << std::endl;
    main_task_factory_->push_back( new core::pack::task::operation::ReadResfile );
  }
	//core::scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
  //energymethodoptions.hbond_options()->decompose_bb_hb_into_pair_energies(true);
  //scorefxn_->set_energy_method_options( energymethodoptions );
}

//setup the packer task and the movemap
void SymFixbbMover::setup_task_mm(core::pose::Pose & pose, core::pack::task::TaskFactoryOP tf  ){
  using namespace core;
  //make the task
  task_ = tf->create_task_and_apply_taskoperations( pose );
  //make movemap to mirror task
  movemap_->set_jump(false);
  movemap_->set_bb(false);
  for( Size ii = 1; ii<= pose.total_residue(); ++ii){
    if( task_->pack_residue( ii ) )
      movemap_->set_chi(ii, true);
    else
      movemap_->set_chi(ii, false);
  }
}//end setup task and movemap


//apply
void SymFixbbMover::apply( core::pose::Pose & pose ) {
  using namespace protocols::simple_moves::symmetry;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

//   //create a task factory: this will create a new PackerTask for each input pose
//   core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
//   main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
//   if ( option[ packing::resfile ].user() ) {
// 		TR << "Reading a resfile for packing." << std::endl;
//     main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
//   }


  //setup symmetric pose
  SetupForSymmetryMoverOP setup_mover = new SetupForSymmetryMover;
  setup_mover->apply(pose);

	//call to helper function to set up the tasks
  setup_task_mm(pose, main_task_factory_);

  //setup the movers
  SymPackRotamersMoverOP sym_pack = new SymPackRotamersMover(scorefxn_, task_);

  //This sequence mover will contain packing for sure, and may contain minimization
  protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
  seq_mover->add_mover( sym_pack );

  //minimize sidechains in task if wanted
  if( option[ minimize_sidechains ] ){
    SymMinMoverOP sym_min = new SymMinMover(movemap_,
																						scorefxn_,
																						option[ OptionKeys::run::min_type ].value(),
																						0.01,
																						true /*use_nblist*/ );

		TR << "Minimizing after packing with:  "<< option[ OptionKeys::run::min_type ].value() << std::endl;
		//add minimization
    seq_mover->add_mover(sym_min);
 }

  //apply mover(s)
  seq_mover->apply(pose);

}//end apply


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  option.add( minimize_sidechains, "Do minimization of side chains after rotamer packing").def(false);

  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go( new SymFixbbMover );

  std::cout << "Done! -------------------------------"<< std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} //end main


