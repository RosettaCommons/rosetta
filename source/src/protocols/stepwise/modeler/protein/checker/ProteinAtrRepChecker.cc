// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility>

static basic::Tracer TR( "protocols.stepwise.modeler.protein.checker.ProteinAtrRepChecker" );

using namespace core;

// AMW TODO:
// Refactor this so that we don't do a multiplication every time we check the screen.
// Instead, change the cutoff values accordingly so they account for unweighted vals.

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace checker {

//Constructor
ProteinAtrRepChecker::ProteinAtrRepChecker( pose::Pose const & pose,
	utility::vector1< Size > const & moving_res_list ):
	moving_res_list_(  moving_res_list )
{
	rep_cutoff_ = 10.0;
	atr_cutoff_ = -0.4; // allows a single H-bond to count.
	initialize_scorefxn();
	get_base_atr_rep_score( pose );
}

//Destructor
ProteinAtrRepChecker::~ProteinAtrRepChecker() = default;

///////////////////////////////////////////
void
ProteinAtrRepChecker::initialize_scorefxn(){
	// Bare minimum to check for contact (fa_atr) but not clash (fa_rep)
	atr_rep_screening_scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_atr , 0.23 );
	atr_rep_screening_scorefxn_->set_weight( core::scoring::fa_rep , 0.12 );
}

///////////////////////////////////////////
void
ProteinAtrRepChecker::get_base_atr_rep_score( core::pose::Pose const & pose ){

	using namespace core::scoring;
	using namespace core::pose;

	pose::Pose base_pose_screen = pose; //hard copy
	base_pose_screen.remove_constraints(); // floating point errors if coordinate constraints are in there.
	split_pose( base_pose_screen, moving_res_list_ );

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// if enabling modeler of loop takeoff, should recognize this and virtualize O at takeoff and H at landing.
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	( *atr_rep_screening_scorefxn_ )( base_pose_screen );
	EnergyMap const & energy_map = base_pose_screen.energies().total_energies();
	base_atr_score_ = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[ scoring::fa_atr ]; //
	base_rep_score_ = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[ scoring::fa_rep ];
}

///////////////////////////////////////////////////////////////////////////////
bool
ProteinAtrRepChecker::check_screen( pose::Pose & current_pose_screen ){
	using namespace core::scoring;
	using namespace ObjexxFCL;

	( *atr_rep_screening_scorefxn_ )( current_pose_screen );

	EnergyMap const & energy_map = current_pose_screen.energies().total_energies();
	Real rep_score = atr_rep_screening_scorefxn_->get_weight( fa_rep ) * energy_map[scoring::fa_rep];
	Real atr_score = atr_rep_screening_scorefxn_->get_weight( fa_atr ) * energy_map[scoring::fa_atr];

	delta_rep_score_ = rep_score - base_rep_score_;
	delta_atr_score_ = atr_score - base_atr_score_;

	// if ( delta_rep_score_ < (  -0.1 ) ){
	//  std::string const message = "delta_rep_score_ = " + string_of( delta_rep_score_ ) + " rep_score = " + string_of( rep_score ) + " base_rep_score = " + string_of( base_rep_score_ );
	//  current_pose_screen.dump_pdb( "PROBLEM.pdb" );
	//  utility_exit_with_message( "delta_rep_score_ < (  -0.1 ), " + message );
	// }

	// if ( delta_atr_score_ > (  +0.1 ) ){
	//  std::string const message = "delta_atr_score_ = " + string_of( delta_atr_score_ ) + " atr_score = " + string_of( atr_score ) + " base_atr_score = " + string_of( base_atr_score_ );
	//  current_pose_screen.dump_pdb( "PROBLEM.pdb" );
	//  utility_exit_with_message( "delta_atr_score_ > (  +0.1 ), " + message );
	// }

	Real actual_rep_cutoff = rep_cutoff_; //default
	bool pass_rep_screen = ( delta_rep_score_ < actual_rep_cutoff );
	bool pass_atr_screen = ( delta_atr_score_ < atr_cutoff_ );

	bool pass_atr_rep_screen = pass_rep_screen && pass_atr_screen;
	return pass_atr_rep_screen;
}


} //checker
} //protein
} //modeler
} //stepwise
} //protocols

