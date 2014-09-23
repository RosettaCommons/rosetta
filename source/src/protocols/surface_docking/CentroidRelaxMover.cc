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
/// @author Michael Pacella (mpacella88@gmail.com)

// Unit Headers
#include <protocols/surface_docking/CentroidRelaxMover.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Utility Headers
#include <utility/exit.hh>

//Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.SurfaceDocking.CentroidRelaxMover" );


namespace protocols {
namespace surface_docking {

using namespace core;
using namespace protocols::moves;

//constructor
CentroidRelaxMover::CentroidRelaxMover() : Mover()
	{
		init();
	}
	
CentroidRelaxMover::CentroidRelaxMover(CentroidRelaxMover const & src) : Mover(src)
	{
		copy_data(*this, src);
	}
	
//destructor
CentroidRelaxMover::~CentroidRelaxMover() {}
	
protocols::moves::MoverOP
CentroidRelaxMover::clone() const
	{
		return protocols::moves::MoverOP( new CentroidRelaxMover(*this) );
	}
	
protocols::moves::MoverOP
CentroidRelaxMover::fresh_instance() const
	{
		return protocols::moves::MoverOP( new CentroidRelaxMover() );
	}
	
void CentroidRelaxMover::show(std::ostream & output) const
	{
		using namespace std;
		
		Mover::show(output);  // name, type, tag
	}


void CentroidRelaxMover::setup_defaults()
	{
	TR << "Setting Defaults" << std::endl;
	// setting MoveMaps
	// Setting up default score function -->place to add modified score6 from
	// rosetta++
	score_low_res_ =
			scoring::ScoreFunctionFactory::create_score_function("RS_centroid.wts");//("cen_std.wts");
	small_min_type_="linmin";
	shear_min_type_="dfpmin";
	}


void CentroidRelaxMover::setup_movers( const core::pose::Pose & pose)
	{
	// Setting Common Parameters for movers
	core::kinematics::MoveMapOP move_map( new core::kinematics::MoveMap() );
	move_map->set_bb( true );
	Real tolerance=0.01;// for minimizer

	//Setting up smallTrialMover
	//Creating smallMove; temperature_ : explanation;
	//          nnmoves_ : number of residues to move(Robin thinks its a
	//          good idea to change this based on length of peptide !!!!
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover(move_map ,kT_,nmoves_/2) );
	small_mover->angle_max(30); // max angle deviation.
	//Default Values copied from  MinMover.cc
	simple_moves::MinMoverOP small_min_mover( new simple_moves::MinMover(move_map ,
                      score_low_res_,small_min_type_,tolerance,true,false,false) );
        //smallsequenceMover_ =
	///     new moves::SequenceMover(smallmover_,smallminmover_);
	moves::SequenceMoverOP small_sequence_mover( new moves::SequenceMover() );
	small_sequence_mover->add_mover(small_mover);
	small_sequence_mover->add_mover(small_min_mover);

	//Setting up shearTrialMover
	//Creating smallMove; temperature_ : explanation;
	//nnmoves_ : number of residues to move
	simple_moves::ShearMoverOP shear_mover( new simple_moves::ShearMover(move_map ,kT_,nmoves_) );
	shear_mover->angle_max(30);
	//Default Values copied from  MinMover.cc
	simple_moves::MinMoverOP shear_min_mover( new simple_moves::MinMover(move_map ,score_low_res_,
	                 shear_min_type_,tolerance,true,false,false) );
	
	moves::SequenceMoverOP shear_sequence_mover( new moves::SequenceMover() );
	shear_sequence_mover->add_mover(shear_mover);
	shear_sequence_mover->add_mover(shear_min_mover);

	TR << "Finalizing Movers" << std::endl;
	// To init the final part of min_trial_small & mini_trial_shear mover
	// Creating a global mc object that can be used by both
	// small_min_trial_mover() and shear_min_trial_mover()
	monte_carlo_ = moves::MonteCarloOP( new moves::MonteCarlo(pose,*score_low_res_,kT_) );
	//smallTrialMover
	TR << "Finalizing small Movers" << std::endl;
	small_trial_min_mover_ = moves::TrialMoverOP( new moves::TrialMover(small_sequence_mover ,monte_carlo_) );
	//shearTrialMover
	TR << "Finalizing shear Movers" << std::endl;
	shear_trial_min_mover_ = moves::TrialMoverOP( new moves::TrialMover(shear_sequence_mover ,monte_carlo_) );
	}

void CentroidRelaxMover::apply(pose::Pose & pose)
	{
	setup_movers(pose);
	
	for ( Size outer_loop = 1; outer_loop <= outer_loop_cycles_; ++outer_loop )
	{
			TR << "CentroidRelax OuterLoop Execution Number:"<< outer_loop << std::endl;
			monte_carlo_->recover_low(pose);
			shear_trial_min_mover_->apply(pose);
			for ( Size inner_loop = 1; inner_loop <=inner_loop_cycles_; ++inner_loop )
			{
				TR << "CentroidRelax InnerLoop Execution Number:"<< inner_loop << std::endl;
				small_trial_min_mover_->apply(pose);
			}
	}
	monte_carlo_->show_counters();
	monte_carlo_->recover_low(pose);
	}

void CentroidRelaxMover::set_nmoves( core::Size const nmoves_in )
	{
		nmoves_ = nmoves_in;
	}
	
void CentroidRelaxMover::set_inner_loop_cycles(const core::Size setting)
	{
		inner_loop_cycles_ = setting;
	}

void CentroidRelaxMover::set_outer_loop_cycles(const core::Size setting)
	{
		outer_loop_cycles_ = setting;
	}
std::string CentroidRelaxMover::get_name() const
	{
	return "CentroidRelaxMover";
	}
	
void CentroidRelaxMover::init()
	{
		nmoves_=2;
		kT_=0.5;
		Mover::type( "CentroidRelaxMover");
		TR << "CentroidRelaxMover Constructor Called" << std::endl;
		setup_defaults();
	}

void CentroidRelaxMover::copy_data(CentroidRelaxMover object_to_copy_to, CentroidRelaxMover object_to_copy_from)
	{
		object_to_copy_to.nmoves_ = object_to_copy_from.nmoves_;
		object_to_copy_to.kT_ = object_to_copy_from.kT_;
		object_to_copy_to.score_low_res_ = object_to_copy_from.score_low_res_;
		object_to_copy_to.small_min_type_ = object_to_copy_from.small_min_type_;
		object_to_copy_to.shear_min_type_ = object_to_copy_from.shear_min_type_;
		object_to_copy_to.monte_carlo_ = object_to_copy_from.monte_carlo_;
		object_to_copy_to.small_trial_min_mover_ = object_to_copy_from.small_trial_min_mover_;
		object_to_copy_to.shear_trial_min_mover_ = object_to_copy_from.shear_trial_min_mover_;
		object_to_copy_to.inner_loop_cycles_ = object_to_copy_from.inner_loop_cycles_;
		object_to_copy_to.outer_loop_cycles_ = object_to_copy_from.outer_loop_cycles_;
	}

}	//surfaceDockingProtocol

}	//protocol
