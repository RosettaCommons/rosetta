// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//      under license.
// (c) The Rosetta software is developed by the contributing members of the
//      Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
//      Questions about this can be
// (c) addressed to University of Washington UW TechTransfer,
//                                             email: license@u.washington.edu.

/// @file FullatomRelaxMover.cc
/// @author Robin A Thottungal

//Unit header
#include <protocols/surface_docking/FullatomRelaxMover.hh>

// Package header
#include <protocols/surface_docking/SlideIntoSurface.hh>
#include <protocols/surface_docking/SurfaceOrientMover.hh>

// Project headers

// AUTO-REMOVED #include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockMCMProtocol.hh>
// AUTO-REMOVED #include <protocols/docking/DockFilters.hh> // get error if you did not include
#include <protocols/moves/RigidBodyMover.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/BackboneMover.fwd.hh>
#include <protocols/moves/MinMover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>

// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>

// C++ Headers
#include <utility/exit.hh>
// AUTO-REMOVED #include <basic/prof.hh>

// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/BackboneMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.SurfaceDocking.FullatomRelaxMover");


namespace protocols {
namespace surface_docking {

using namespace core;


//constructor
FullatomRelaxMover::FullatomRelaxMover() : Mover(){
	// initilized here to avoid warning!
	nmoves_=1;
	temperature_=0.5;
	kT_=0.5;
	Mover::type("FullatomRelaxMover");
	TR << "FullatomRelaxMover Constructor Called" << std::endl;
	setup_defaults();
}

//destructor
FullatomRelaxMover::~FullatomRelaxMover() {}


void FullatomRelaxMover::setup_defaults(){
	TR << "Setting Defaults" << std::endl;
	// setting MoveMaps
	//Setting up default scorefunction --> rosetta++ modified score12
	score_high_res_ = scoring::ScoreFunctionFactory::create_score_function
                                                    ( "standard","score12" );
	smallmin_type_="linmin";
	shearmin_type_="dfpmin";
	benchmark_=false;
	encounter_cycle_=2; // need to get replaced by a random number
	setupMovers();
	}

void FullatomRelaxMover::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	benchmark_ = option[ OptionKeys::run::benchmark ]();

}

void FullatomRelaxMover::setupMovers(){
	// Setting Common Parameters for movers
	moveMapOP_=new core::kinematics::MoveMap;
	//moveMapOP_->set_bb( true );
	Real tolerance=0.01;// for minimizer
	//Setting up smallTrialMover
	//Creating smallMove; temperature_ : explanation;
	//          nnmoves_ : number of residues to move(Robin thinks its a
	//          good idea to change this based on length of peptide !!!!
	smallmover_=new moves::SmallMover(moveMapOP_,temperature_,nmoves_);
	smallmover_->angle_max(30); // max angle deviation.
	//Default Values copied from  MinMover.cc
	smallminmover_= new moves::MinMover(moveMapOP_,
			score_high_res_,smallmin_type_,tolerance,true,false,false);
	//smallsequenceMover_ =
	///     new moves::SequenceMover(smallmover_,smallminmover_);
	smallsequenceMover_= new moves::SequenceMover;
	smallsequenceMover_->add_mover(smallmover_);
	smallsequenceMover_->add_mover(smallminmover_);

	//Setting up shearTrialMover
	//Creating smallMove; temperature_ : explanation;
	//									nnmoves_ : number of residues to move
	shearmover_=new moves::ShearMover(moveMapOP_,temperature_,nmoves_/2);
	shearmover_->angle_max(30); // need to enable after testing!!
	//Default Values copied from  MinMover.cc
	shearminmover_= new moves::MinMover(moveMapOP_,score_high_res_,
					shearmin_type_,tolerance,true,false,false);
	//shearsequenceMover_ =
	//     new moves::SequenceMover(shearmover_,shearminmover_);
	shearsequenceMover_= new moves::SequenceMover;
	shearsequenceMover_->add_mover(shearmover_);
	shearsequenceMover_->add_mover(shearminmover_);

}


void FullatomRelaxMover::FinalizeMovers(pose::Pose & pose)
{
	TR << "Finalizing Movers" << std::endl;
	// Setting the move map for the peptide
	// For setting movemap in fullatom
	Size protein_startseqnum=0;
	for (Size i=1; i<=pose.total_residue(); ++i){
		if (pose.residue_type(i).is_protein()){
			protein_startseqnum=i;
			break;
		}

	}
	TR<<"Protein Start Location:"<<protein_startseqnum<<std::endl;
	moveMapOP_->set_bb_true_range(protein_startseqnum,pose.total_residue());

	//core::pose::PDBInfoOP pdbinfo = pose.pdb_info();
	// Hard coded the chainID hear, need to come up with a smart way of
	// getting this info Stupid way !! temp fix @Robin
	//TR<<"Movemap Range Begins here: "<<pdbinfo->pdb2pose('B',1)<<std::endl;
	//moveMapOP_->set_bb_true_range(pdbinfo->pdb2pose('B',1),pose.total_residue());
	// Setting the jump between the protein and the surface movable
	moveMapOP_->set_jump(pose.num_jump(),true);

	// To init the final part of the min_trial_small & mini_trial_shear mover
	// Creating a global mc object that can be used by both
	// small_min_trial_mover() and shear_min_trial_mover()
	monteCarlo_ = new moves::MonteCarlo(pose,*score_high_res_,kT_);
	//smallTrialMover
	smallmonteCarlo_ = new moves::MonteCarlo(pose,*score_high_res_,kT_);
	small_trial_min_mover_ =
			new moves::TrialMover(smallsequenceMover_,monteCarlo_);
	//shearTrialMover
	shearmonteCarlo_ = new moves::MonteCarlo(pose,*score_high_res_,kT_);
	shear_trial_min_mover_ =
			new moves::TrialMover(shearsequenceMover_,monteCarlo_);
}

void FullatomRelaxMover::set_smallmovesize(Size scale){
	smallmover_->angle_max(30/scale);
}

void FullatomRelaxMover::set_ljrepulsion_weight(Real weight_scale){
	score_high_res_->set_weight(core::scoring::fa_rep,weight_scale);
}

void FullatomRelaxMover::set_ecounter(Size ecount){
	encounter_=ecount;
}


void FullatomRelaxMover::apply(pose::Pose & pose){
using namespace docking;
using namespace moves;
FinalizeMovers(pose);
//AddPyMolObserver(pose, false);
TR << "Starting FullatomRelax" << std::endl;
core::Size upper_limit=3;
// object for slide into contact
FaSlideIntoSurface surfaceContact( pose.num_jump());
for ( Size j = 1; j <=upper_limit; ++j ){ // Original is 5
    if ( encounter_ > encounter_cycle_ ){
		//slide into contact
		surfaceContact.apply( pose );
		}
    smallminmover_->min_type("dfpmin");
    shearminmover_->min_type("dfpmin");
    small_trial_min_mover_->apply(pose); // dfpmin
    shear_trial_min_mover_->apply(pose);  //dfpmin
    //crank mover is missing
    for ( Size k = 1; k <= upper_limit; ++k ){//k<=5 original value;k=2 test
        smallminmover_->min_type("linmin");
        shearminmover_->min_type("linmin");
        small_trial_min_mover_->apply(pose); //linmin
        shear_trial_min_mover_->apply(pose);  //linmin
        //crank mover is missing
        }
    if ( encounter_ == encounter_cycle_ && j%upper_limit == 0 ){//original 5
        //Write the solutionState Structure
        // Random Orient the Partner (make sure this is the peptide)
        TR<<"RigidBodyRandomizeMover"<<std::endl;
        Size rb_jump_=pose.num_jump(); //default value
        RigidBodyRandomizeMover rmover( pose, rb_jump_, partner_upstream );
        rmover.apply( pose );
        //Axis Spin
        RigidBodySpinMover smover( rb_jump_ );
        smover.apply( pose );
        // SurfaceOrient Mover
        surface_docking::SurfaceOrientMoverOP sf=
                                   new surface_docking::SurfaceOrientMover();
        sf->apply(pose);
        surfaceContact.apply( pose );
        monteCarlo_->reset(pose);
        //EncounterComplexFormed = true;
        }
		if ( encounter_ >= encounter_cycle_ && j%upper_limit== 0 ){
                                                                   //original 5
        TR<<"Started HighResolution Surface Docking "<<std::endl;
        moveMapOP_->set_chi(true);
        DockingHighResOP dockmcm =new DockMCMProtocol(pose.num_jump());
        dockmcm->apply(pose);
        //minimize_set_vary_chi(true); --> need to identify
        //docking_monte_carlo_minimize (6,"dfpmin",score12,0.1,rotmag,15,1.0);
        }
	}


}

void FullatomRelaxMover::set_nmoves( core::Size const nmoves_in ){
	nmoves_ = nmoves_in;
    }

std::string FullatomRelaxMover::get_name() const {
    return "FullatomRelaxMover";
   }



}//surfaceDockingProtocol

} //protocol
