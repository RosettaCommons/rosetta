// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
//     vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file CentroidRelaxMover.cc
/// @author Robin A Thottungal  (rathottungal@gmail.com)

// Unit Headers
#include <protocols/surface_docking/CentroidRelaxMover.hh>

// Package Headers

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/jobdist/Jobs.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

//Utility Headers
#include <utility/exit.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.SurfaceDocking.CentroidRelaxMover");


namespace protocols {
namespace surface_docking {

using namespace core;
using namespace protocols::moves;

//constructor
CentroidRelaxMover::CentroidRelaxMover() : Mover(){
	// initilized here to avoid warning!
	nmoves_=2;
	temperature_=0.5;
	kT_=0.5;
	Mover::type( "CentroidRelaxMover");
	TR << "CentroidRelaxMover Constructor Called" << std::endl;
	setup_defaults();
	}

//destructor
CentroidRelaxMover::~CentroidRelaxMover() {}

void CentroidRelaxMover::setup_defaults(){
	TR << "Setting Defaults" << std::endl;
	// setting MoveMaps
	// Setting up default score function -->place to add modified score6 from
	// rosetta++
	score_low_res_ =
			scoring::ScoreFunctionFactory::
					create_score_function("RS_centroid.wts");//("cen_std.wts");
	smallmin_type_="linmin";
	shearmin_type_="dfpmin";
	benchmark_=false;
	setupMovers();
	init_from_options();
	}

void CentroidRelaxMover::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	benchmark_ = option[ OptionKeys::run::benchmark ]();

}

void CentroidRelaxMover::setupMovers(){
	// Setting Common Parameters for movers
	moveMapOP_=new core::kinematics::MoveMap;
	moveMapOP_->set_bb( true );
	Real tolerance=0.01;// for minimizer

	//Setting up smallTrialMover
	//Creating smallMove; temperature_ : explanation;
	//          nnmoves_ : number of residues to move(Robin thinks its a
	//          good idea to change this based on length of peptide !!!!
	smallmover_=new simple_moves::SmallMover(moveMapOP_,temperature_,nmoves_/2);
	smallmover_->angle_max(30); // max angle deviation.
	//Default Values copied from  MinMover.cc
	smallminmover_= new simple_moves::MinMover(moveMapOP_,
                      score_low_res_,smallmin_type_,tolerance,true,false,false);
        //smallsequenceMover_ =
	///     new moves::SequenceMover(smallmover_,smallminmover_);
	smallsequenceMover_= new moves::SequenceMover;
	smallsequenceMover_->add_mover(smallmover_);
	smallsequenceMover_->add_mover(smallminmover_);

	//Setting up shearTrialMover
	//Creating smallMove; temperature_ : explanation;
	//nnmoves_ : number of residues to move
	shearmover_=new simple_moves::ShearMover(moveMapOP_,temperature_,nmoves_);
	shearmover_->angle_max(30);
	//Default Values copied from  MinMover.cc
	shearminmover_= new simple_moves::MinMover(moveMapOP_,score_low_res_,
	                 shearmin_type_,tolerance,true,false,false);
	//shearsequenceMover_ =
        //new moves::SequenceMover(shearmover_,shearminmover_);
	shearsequenceMover_= new moves::SequenceMover;
	shearsequenceMover_->add_mover(shearmover_);
    shearsequenceMover_->add_mover(shearminmover_);

}


void CentroidRelaxMover::FinalizeMovers(pose::Pose & pose){
	TR << "Finalizing Movers" << std::endl;
	// To init the final part of min_trial_small & mini_trial_shear mover
	// Creating a global mc object that can be used by both
	// small_min_trial_mover() and shear_min_trial_mover()
	monteCarlo_ = new moves::MonteCarlo(pose,*score_low_res_,kT_);
	//smallTrialMover
	TR << "Finalizing small Movers" << std::endl;
	smallmonteCarlo_ = new moves::MonteCarlo(pose,*score_low_res_,kT_);
	small_trial_min_mover_ =
		   new moves::TrialMover(smallsequenceMover_,monteCarlo_);
	//shearTrialMover
	TR << "Finalizing shear Movers" << std::endl;
	shearmonteCarlo_ = new moves::MonteCarlo(pose,*score_low_res_,kT_);
	shear_trial_min_mover_ =
	           new moves::TrialMover(shearsequenceMover_,monteCarlo_);
	}

void CentroidRelaxMover::apply(pose::Pose & pose){
	FinalizeMovers(pose);
	TR << "Starting CentroidRelax" << std::endl;
	Size innerLoopCycle= pose.total_residue();
	Size outerLoopCycle=6;
	for ( Size outerLoop = 1; outerLoop <= outerLoopCycle; ++outerLoop ){
            TR << "OuterLoop Execution Number:"<< outerLoop << std::endl;
            monteCarlo_->recover_low(pose);
            shear_trial_min_mover_->apply(pose);
            for ( Size innerLoop = 1; innerLoop <=innerLoopCycle; ++innerLoop ){
                TR << "InnerLoop Execution Number:"<< innerLoop << std::endl;
                small_trial_min_mover_->apply(pose);
            }
	}
	monteCarlo_->show_counters();
	monteCarlo_->recover_low(pose);
}

void CentroidRelaxMover::set_nmoves( core::Size const nmoves_in ){
	nmoves_ = nmoves_in;
	}

std::string CentroidRelaxMover::get_name() const {
	return "CentroidRelaxMover";
	}



}	//surfaceDockingProtocol

}	//protocol
