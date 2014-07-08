// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/stranges/peptide_minimizer.cc
/// @brief does stuff
/// @author Ben Stranges

// Unit headers
#include <devel/init.hh>

//  Core library
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>


//  Movers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh> //Small Mover
#include <protocols/moves/MoverContainer.hh> //SequenceMover
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/docking/DockingProtocol.hh>

// Job distributor
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
//#include <core/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>



// C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>


//using core::util::Error;
//using core::util::Warning;
static basic::Tracer TR("PepMin");
static basic::Tracer TRdebug("PepMin.debug");

using namespace core;
using namespace std;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options;
using namespace basic::options::OptionKeys;

// application specific options

RealOptionKey const bbhbond_weight_mod( "bbhbond_weight_mod" ); //amount to upweight hbond E
IntegerOptionKey const ncycles( "ncycles" ); // number of minimization cycles


// mover deffinition
class PepMinMover : public Mover {
public:

  PepMinMover();

  virtual void apply( core::pose::Pose& pose );

  virtual MoverOP clone() const {
    return new PepMinMover( *this );
  }

  virtual MoverOP fresh_instance() const {
    return clone();
  }

	virtual std::string get_name() const{
		return "PepMinMover";
	}

  virtual void print_energy(std::string const &  pose_name,
			    core::Real & e_value,
			    core::Size & bb_hbonds );

private:

  core::Real hbond_multiplier_;
  core::scoring::ScoreFunctionOP scorefxn_;
  int ncycles_;


};
//constructor
PepMinMover::PepMinMover() {
  hbond_multiplier_ = option[ bbhbond_weight_mod ];
  scorefxn_ = scoring::get_score_function();
  ncycles_ = option[ncycles];

  //set up for scoring
  TR<< "Up/down weighting lr_bb by: " << hbond_multiplier_ << std::endl;
  core::scoring::methods::EnergyMethodOptions emoptions(scorefxn_->energy_method_options());
  emoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
  scorefxn_->set_energy_method_options( emoptions );
  //upweight bb_hbonding
  //Real sr_wt( (*scorefxn_).get_weight(core::scoring::hbond_sr_bb) );
  Real lr_wt( (*scorefxn_).get_weight(core::scoring::hbond_lr_bb) );
  //(*scorefxn_).set_weight(core::scoring::hbond_sr_bb, ( sr_wt * hbond_multiplier_ ) );
  (*scorefxn_).set_weight(core::scoring::hbond_lr_bb, ( lr_wt * hbond_multiplier_ ) );
}


///////////////////////////////////////////////////////////////
void
PepMinMover::print_energy(
	     std::string const &  pose_name,
	     core::Real & e_value,
	     core::Size & bb_hbonds )
{
  using namespace std;

  cout << "FILE:  " << setw(25) << pose_name << "   "
       << "HBONDS:  " << setw(3) << bb_hbonds << "    "
       << "TOTAL_E:  " << setw(6) << setprecision(3) << e_value << "    "
       << "E/HBOND:  " << setw(6) << setprecision(3) << (e_value / bb_hbonds ) <<endl;
}
///////////////////////////////////////////////////////////////

//begin mover apply
void PepMinMover::apply (pose::Pose& pose ) {


  Real const mc_temp( 0.8 );

  //Job stuff
  protocols::jd2::JobOP const job_me( jd2::JobDistributor::get_instance()->current_job() );
  std::string job_name (jd2::JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

  //utility::file::FileName posename (pose.pdb_info()->name());

	//output the foldtree
	TRdebug << "FoldTree in apply. \n" << pose.fold_tree()<< std::endl;

  //make a list of hydrogen bonds
  core::scoring::hbonds::HBondSet hbond_set;

  //make a movemap
  core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  //set move map conditions
  movemap->set_jump(true);
  Size chainAend = pose.conformation().chain_end(1);
	TR << "chain A end: " << chainAend << endl;
  Size nres= pose.total_residue();
  //set bb to moveable for peptide
  for (Size i = chainAend + 1; i <= nres; ++i){
    movemap->set_bb(i, true);
    TRdebug << "bb residue to optimize: " << i << std::endl;
  }

//   //set up for scoring
//   TR<< "Up/down weighting lr_bb by: " << hbond_multiplier_ << std::endl;
//   core::scoring::methods::EnergyMethodOptions emoptions(scorefxn_->energy_method_options());
//   emoptions.decompose_bb_hb_into_pair_energies(true);
//   scorefxn_->set_energy_method_options( emoptions );
//   //upweight bb_hbonding
//   //Real sr_wt( (*scorefxn_).get_weight(core::scoring::hbond_sr_bb) );
//   Real lr_wt( (*scorefxn_).get_weight(core::scoring::hbond_lr_bb) );
//   //(*scorefxn_).set_weight(core::scoring::hbond_sr_bb, ( sr_wt * hbond_multiplier_ ) );
//   (*scorefxn_).set_weight(core::scoring::hbond_lr_bb, ( lr_wt * hbond_multiplier_ ) );


  //make the movers
  //docking
  Size jump (1);
  pose::PoseOP pose_op ( new pose::Pose( pose ));
  docking::DockingProtocolOP dock_mover= new docking::DockingProtocol ( jump, false /*full_atom*/, true /*local refine*/ );
  dock_mover->set_highres_scorefxn( scorefxn_ );
  dock_mover->moves::Mover::set_input_pose( pose_op );
  dock_mover->moves::Mover::set_native_pose( pose_op );

  //minimization
  protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );
  TR<< "Minimizing with: " << option[ OptionKeys::run::min_type ].value() << std::endl;
  //backbone moves
  Size nmoves (5);
  protocols::simple_moves::SmallMoverOP small_mover = new protocols::simple_moves::SmallMover(movemap, mc_temp, nmoves);
  small_mover->set_angles(10.0);
  protocols::simple_moves::ShearMoverOP shear_mover = new protocols::simple_moves::ShearMover(movemap, mc_temp, nmoves);
  shear_mover->set_angles(10.0);
  RandomMoverOP bb_random_mover( new moves::RandomMover() );
  bb_random_mover->add_mover( small_mover );
  bb_random_mover->add_mover( shear_mover );
  //monte carlo mover
  MonteCarloOP mc_mover = new MonteCarlo ( pose, *scorefxn_, mc_temp );
  //make Sequence mover that makes small move then  minimizes
  SequenceMoverOP bb_min_mover = new SequenceMover;
  bb_min_mover->add_mover( bb_random_mover );
  bb_min_mover->add_mover( min_mover );
  //mc trials (what is acutally applied)
  TrialMoverOP trial_mover = new TrialMover(bb_min_mover, mc_mover);


  /////////////////////////////////////////////////////
  //apply movers!!!
  for(Size i =1; i <= 5; ++i){  //randomize peptide a bit
    bb_random_mover->apply( pose );
  }
  dock_mover->apply( pose );  //docks  before searching space
	//Allow this mover to see the status of DockingProtocol because it can fail and that failure needs to be known by this mover.
	TR << "DOCKING MOVE STATUS: (0 = MS_SUCESS) " << dock_mover->get_last_move_status() << endl;
	set_last_move_status( dock_mover->get_last_move_status() );

	//cout << "MOVE STATUS: " << get_last_move_status() << endl;
	if(get_last_move_status() != protocols::moves::MS_SUCCESS){
		TR << "Docking Failed. Restarting." << std::endl;
		return;
	}
  TR << "Searching space..." << std::endl;
  Size const nsteps (ncycles_);
  for (Size ii = 1; ii <= nsteps; ++ii){
    trial_mover->apply( pose );
    if ( ii % (nsteps / 10) == 0)
      TR << job_name << " peptide step " << ii  << " of " << nsteps << "  score:  "<<(*scorefxn_)(pose)<<std::endl;
  }
  //recover best structure
  TR << "Output lowest energy backbone of peptide..." << std::endl;
  mc_mover->recover_low( pose );
///////////////////////////////////////////////

  //figure out energy statistics
  Size n_hbonds (0);
  Real hbond_E (0.0);
  Real hb_wt( (*scorefxn_).get_weight(core::scoring::hbond_lr_bb) );
  //find backbone hbonds for pdb
  hbond_set.clear(); //probably don't need anymore
  pose.update_residue_neighbors();
  core::scoring::hbonds::fill_hbond_set(
					pose,
					false /*calc_deriv*/,
					hbond_set,
					true /*bb only*/ );
  //call to try to resize bb_don/accept arrays
  //need this for everything to work right
  hbond_set.setup_for_residue_pair_energies(pose);
  for (Size i = chainAend + 1; i <= nres; ++i){
    core::scoring::EnergyMap emap( pose.energies().residue_total_energies( i ));
    //how many hbonds?
    if (hbond_set.acc_bbg_in_bb_bb_hbond(i) == true)
      ++n_hbonds;
    if (hbond_set.don_bbg_in_bb_bb_hbond(i) == true)
      ++n_hbonds;
    hbond_E += (hb_wt * emap[core::scoring::hbond_lr_bb]);
  }
  print_energy( job_name,
		hbond_E,
		n_hbonds);
  //pass info to JD2
  job_me->add_string_real_pair("n_hbonds", n_hbonds );
  job_me->add_string_real_pair( "hbondE", hbond_E );
  job_me->add_string_real_pair( "E/Hbond", (hbond_E/n_hbonds) );

}//end apply


// run protocol
int
main( int argc, char * argv [] )
{

	try {

  using namespace protocols;
  using namespace protocols::jd2;
  using namespace protocols::moves;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  option.add( bbhbond_weight_mod, "multiplier for lr bb-bb hbond .").def(1.0);
  option.add( ncycles, "number of perturb/min cycles.").def(100);

  // initialize core
  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go( new PepMinMover );

  std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
