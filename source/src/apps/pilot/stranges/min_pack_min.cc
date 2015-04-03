// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/stranges/min_pack_min.cc
/// @brief does cycles of Minimiziation, Repacking, Minimization for interface analysis later
/// no design by default
/// @author Ben Stranges


// Unit headers
#include <devel/init.hh>

// Project Headers
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// #ifndef NDEBUG
//   #include <core/kinematics/FoldTree.hh>
// #endif

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Job distributor
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>


// Utility Headers
#include <basic/Tracer.hh>


// C++ headers
#include <sstream>
#include <iostream>
#include <string>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "apps.pilot.stranges.min_pack_min" );

using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::jd2;

// application specific options
namespace MinPackMin {
	IntegerVectorOptionKey const rbjump( "rbjump" ); //What jumps to minimize
	BooleanOptionKey const pack_first( "pack_first" ); //pack before min, good if big clashes!
	BooleanOptionKey const no_rbmin( "no_rbmin" ); //do now allow RB min at all
	BooleanOptionKey const min_all_jumps( "min_all_jumps" ); //minimize all jumps
	//BooleanOptionKey const pymolreport( "pymolreport" ); //use pymol mover

}

// mover deffinition
class MinPackMinMover : public Mover {
public:

	MinPackMinMover();

	virtual void apply( core::pose::Pose& pose );

	virtual MoverOP clone() const {
		return MoverOP( new MinPackMinMover( *this ) );
	}

	virtual
	std::string
	get_name() const {
		return "MinPackMinMover";
	}

	virtual	MoverOP	fresh_instance() const {
		return MoverOP( new MinPackMinMover );
	}

private:
	//vector1<Size> default_jump_;
	vector1<int> jump_num_ ;
	bool no_rbmin_, pack_first_;
	bool min_all_;
	scoring::ScoreFunctionOP scorefxn_;
};

//constructor
MinPackMinMover::MinPackMinMover(){
	vector1<Size> default_jump;
	default_jump.push_back(1);
	jump_num_ = option[ MinPackMin::rbjump ].def(default_jump);
	min_all_ = option[ MinPackMin::min_all_jumps ].def(true);
	no_rbmin_ = option[ MinPackMin::no_rbmin].def(false);
	pack_first_ = option[ MinPackMin::pack_first].def(false);
	//pymolreport_ = option[ MinPackMin::pymolreport ];
	scorefxn_ = scoring::get_score_function();
}

//begin mover apply
void MinPackMinMover::apply (pose::Pose& pose ) {

// 	//for pymol viewing
// 	if( pymolreport_ ){
// 		TR << "Adding pymol reporting to pose." << std::endl;
// 		protocols::moves::PyMolObserverOP pymol_ob = 	protocols::moves::AddPyMolObserver(pose, false, 0.33);
// 	}


  //////////////////////////////
  // Handle inputs
  //////////////////////////////

	if(no_rbmin_){
		TR << "Not allowing any rigid body minimization" << std::endl;
	}
	else if (min_all_){
		TR << "Minimizing all jumps in pose." << std::endl;}
	else {
		TR << "Allowing certian rigid body minimization: \n";
		TR << "Rigid body jump defined at number(s): ";
		for(Size ii=1; ii<= jump_num_.size(); ++ii){
			TR << jump_num_[ii] << ", " ;
		}
		TR<< std::endl;
	}

	// get the job info
  protocols::jd2::JobOP const job_me( JobDistributor::get_instance()->current_job() );
	std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

  //define a native pose if given one for RMSD calc later...
  pose::Pose native_pose;
  if (basic::options::option[ in::file::native ].user()){
    core::import_pose::pose_from_pdb( native_pose, basic::options::option[ in::file::native ]());
  }

	// setup for bb-bb pair energies
  scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
  scorefxn_->set_energy_method_options( energymethodoptions );

  // score the input pose for kicks and giggles
  (*scorefxn_)( pose );


// #ifndef NDEBUG
//   TR<< "FoldTree: \n" << pose.fold_tree() << std::endl;
// #endif


  //////////////////////////////////////////////////////
  // Set up packer tasks and move map
  //////////////////////////////////////////////////////

  TaskFactoryOP tf( new TaskFactory() );
  //set common taskfactory options
  tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) );
  tf->push_back( TaskOperationCOP( new operation::IncludeCurrent ) );
	// read a resfile. note that design is impossible
	if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ){
		tf->push_back( TaskOperationCOP( new operation::ReadResfile() ) );
		TR << "Resfile read, note that this application DOES NOT DESIGN!" << std::endl;
	}

  //want to just repack the wt pose, NO design
  RestrictResidueToRepackingOP repack_op( new RestrictResidueToRepacking() );
  for( Size ii = 1; ii<= pose.n_residue(); ++ii) {
    repack_op->include_residue( ii );
  }
  //fill task factory with these restrictions
  tf->push_back( repack_op );

#ifndef NDEBUG
  TR<< "Packer Task: " << *(tf->create_task_and_apply_taskoperations(pose));
#endif

  //set move map to allow bb and sc minimization
  kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
  movemap->set_bb(true);
  movemap->set_chi(true);
	movemap->set_jump( false );
	//setup another movemap that does not contain any RB jumps,
	//use this movemap at the begining  makes no sense to allow RB min before packing
	kinematics::MoveMapOP movemap_initial ( movemap->clone() );
	//only need to set RB jump minimization if there is a jump...
	if( pose.conformation().num_chains() > 1){
		//what jumps will we be minimizing today
		if( no_rbmin_ ){ //no rb min
			movemap->set_jump( false );
			TR << "Setting all jumps in movemap to false!" << std::endl;
		}
		else if(min_all_)
			movemap->set_jump(true);
		else{
			//allow RB minimization for the defined jump in the complex only
			// under option control for mutations that may push the complex away in minimization
			movemap->set_jump(false);
			for(Size ii = 1; ii <= jump_num_.size(); ++ii ){
				if( jump_num_[ii] > (int)pose.num_jump() )
					utility_exit_with_message( "Input declared a jump with too high a number, exiting... \n" );
				else
					movemap->set_jump( jump_num_[ii], true);
			}
		}
	}
  /////////////////////////////////////////////////////////
  // Make MinMovers and PackRotsMover
  ////////////////////////////////////////////////////////

  //first do minimization mover making
	//leave out movemap in initialization
  protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap_initial, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ ) );

  // Make PackRots movers
  protocols::simple_moves::PackRotamersMoverOP repacker( new protocols::simple_moves::PackRotamersMover() );
  repacker->task_factory( tf );
  repacker->score_function( scorefxn_ );


  ///////////////////////////////////////////////////////////////////////////////
  //  Now do actual applys  (wahoo!!!)
  ///////////////////////////////////////////////////////////////////////////////

  //Print out initial scores
  TR << "Initial scores for input: " << job_name  <<": "<< (*scorefxn_)(pose) << std::endl;

	if(pack_first_){
		//Step 2 packing incase minimization not enough in some cases (big clashes)
		TR << "Packing first for: "<< job_name << std::endl;
		repacker->apply(pose);
		//Print out scores after packing
		TR << "Score after first packing " << job_name << ": " <<(*scorefxn_)(pose)  << std::endl;

	}

  // Step 1 minimization... no RB min allowed
	min_mover->movemap(movemap_initial);
  min_mover->apply(pose);

  //Print out scores after minimization
  TR << "Score after first min: " << job_name  <<": "<< (*scorefxn_)(pose) << std::endl;

  //Step 2 packing incase minimization not enough in some cases (big clashes)
  TR << "Packing results for: "<< job_name << std::endl;
  repacker->apply(pose);

  //Print out scores after packing
  TR << "Best score after Packing " << job_name << ": " <<(*scorefxn_)(pose)  << std::endl;

  //Step 3
  // One more round of minimization and output.  Allow RB min if given through options
	min_mover->movemap(movemap);
	min_mover->apply( pose );

	//figure out an rmsd if the native structure is given
	core::Real rms (0.0);
	if (basic::options::option[ in::file::native ].user()){
		// allow superposition because RB min is allowed
		rms = scoring::rmsd_with_super( native_pose, pose , scoring::is_protein_CA ) ;
		job_me->add_string_real_pair("native_rms", rms );
	}
	else{
		rms = scoring::rmsd_with_super( *(job_me->get_pose()), pose, core::scoring::is_protein_CA );
		job_me->add_string_real_pair("rms", rms );
	}

	TR << "Score after final minimization: "  << job_name <<": "<< (*scorefxn_)(pose) << "  rmsd: " << rms << std::endl;

  TR << "Complete..." << std::endl;


}//end apply

// run protocol
int
main( int argc, char * argv [] )
{

	try {

	//using namespace protocols;
	using namespace protocols::jd2;
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

  option.add( MinPackMin::rbjump, "Defines the rigid body jump to split the interface apart by.");
  option.add( MinPackMin::pack_first, "Add packing before initial minimization.");
  option.add( MinPackMin::no_rbmin, "Do not allow rigid body minimization." );
	option.add( MinPackMin::min_all_jumps, "Rigid body minimizationon ALL jumps." );
	//option.add( MinPackMin::pymolreport, "Report to pymol observer").def(false);

	// initialize core
	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new MinPackMinMover ) );

	std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


