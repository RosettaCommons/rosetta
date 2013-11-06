// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/DualMonteCarlo.hh
/// @brief "dual" MonteCarlo wraps MonteCarlo to allow for centroid MonteCarlo scoring of a fullatom pose (or other similar situations)
/// @author Steven Lewis

// Unit Headers
#include <protocols/moves/DualMonteCarlo.hh>
#include <protocols/moves/MonteCarlo.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// ObjexxFCL Headers

// Numeric Headers
//#include <numeric/random/random.hh>  //NO, we are not doing any boltzmann here

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//#include <basic/prof.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.DualMonteCarlo");

namespace protocols {
namespace moves {

/// @details ctor; will call MC ctor too
DualMonteCarlo::DualMonteCarlo(
															 core::pose::Pose const & DMC_pose,
															 core::pose::Pose const & MC_pose,
															 core::scoring::ScoreFunction const & DMC_scorefunction,
															 core::scoring::ScoreFunction const & MC_scorefunction,
															 core::Real const temperature
):
	MC_(MC_pose, MC_scorefunction, temperature)
{
	last_accepted_pose_ = new core::pose::Pose();
	lowest_score_pose_ = new core::pose::Pose();
	score_function_ = DMC_scorefunction.clone();
	reset( DMC_pose, MC_pose );
}

DualMonteCarlo::~DualMonteCarlo() {}

void
DualMonteCarlo::reset( core::pose::Pose const & DMC_pose, core::pose::Pose const & MC_pose )
{
	reset_DMC(DMC_pose);
	MC_.reset(MC_pose);
}

void
DualMonteCarlo::reset_DMC( core::pose::Pose const & DMC_pose ){

	//PROF_START( basic::MC_ACCEPT );
	*last_accepted_pose_ = DMC_pose;
	//PROF_STOP( basic::MC_ACCEPT );

	(*score_function_)( *last_accepted_pose_ );
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );

	//PROF_START( basic::MC_ACCEPT );
	*lowest_score_pose_ = *last_accepted_pose_;
	//PROF_STOP( basic::MC_ACCEPT );
}


/// @brief Copies the lowest_score_ pose into the input pose, as well as into
/// last_accepted_pose_.  Copies action for underlying MC but does not return that pose.
void
DualMonteCarlo::recover_low( core::pose::Pose & DMC_pose ){
	core::pose::Pose toss;
	recover_low(DMC_pose, toss); //throw away the other one
}

/// @brief Copies the lowest_score_ pose into the input pose, as well as into
/// last_accepted_pose_.  Always copies action for underlying MC.  This version returns the MC's pose.
void
DualMonteCarlo::recover_low( core::pose::Pose & DMC_pose, core::pose::Pose & MC_pose ){
	( DMC_pose ) = ( *lowest_score_pose_ );
	*last_accepted_pose_ = *lowest_score_pose_ ;
	MC_.recover_low( MC_pose );
}

/////////////////////////////////////////////////////////////////////////////
void
DualMonteCarlo::show_scores() const
{
	TR << "last_accepted_score, lowest_score: " << last_accepted_score() << ", " << lowest_score() << std::endl;
}

// /// @brief change the scorefxn,  re-scores last-accepted and lowest-score pose.  Will NOT change their identities.
// void
// DualMonteCarlo::score_function( core::scoring::ScoreFunction const & DMC_scorefunction )
// {
// 	score_function_ = new ScoreFunction(DMC_scorefunction);

// 	(*score_function_)( *lowest_score_pose_ );
// 	score_function_->accumulate_residue_total_energies( *lowest_score_pose_ );

// 	(*score_function_)( *last_accepted_pose_ );
// 	score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );
// }

// /// @brief change the scorefxns,  re-scores last-accepted and lowest-score pose, calls MC's score_function, and may change the identity of the poses if MC does too.
// void
// DualMonteCarlo::score_function(
// 															 core::scoring::ScoreFunction const & DMC_scorefunction,
// 															 core::scoring::ScoreFunction const & MC_scorefunction)


/////////////////////////////////////////////////////////////////////////////
//////////////////////////////
//
// return true for an accept, false otherwise
//
// 	mc_accepted
// 		3 = accepted:score beat low score and last_accepted score
// 		2 = accepted:score beat last_accepted score
// 		1 = thermally accepted: score worse than last_accepted score
// 		0 = not accepted
/// @details remember that this function DOES NOT CARE what scores the DMC_poses have!
bool
DualMonteCarlo::boltzmann(
													core::pose::Pose & DMC_pose,
													core::pose::Pose & MC_pose,
													std::string const & move_type// = "unk"
)
{
	const bool accepted(MC_.boltzmann(MC_pose, move_type));

	switch (MC_.mc_accepted()) {
	case MCA_accepted_score_beat_low:
		reset_DMC(DMC_pose);
		break;
	case MCA_accepted_score_beat_last:
	case MCA_accepted_thermally:
		//PROF_START( basic::MC_ACCEPT );
		*last_accepted_pose_ = DMC_pose;
		//PROF_STOP( basic::MC_ACCEPT );
		(*score_function_)( *last_accepted_pose_ );
		/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );
		break;
	case MCA_rejected:
		DMC_pose = *last_accepted_pose_;
		break;
	default:
		break;
	}
	return accepted;
}


// // for recovering from a checkpoint
// void
// DualMonteCarlo::set_last_accepted_pose( Pose const & pose )
// {
// 	*last_accepted_pose_ = pose;
// }

// // for recovering from a checkpoint
// void
// DualMonteCarlo::set_lowest_score_pose( Pose const & pose )
// {
// 	*lowest_score_pose_ = pose;
// }

// void
// DualMonteCarlo::attach_observer_to_last_accepted_conformation(
// 	core::conformation::ConformationObserverOP observer
// )
// {
// 	last_accepted_pose_->conformation().attach_observer( observer );
// }

// void
// DualMonteCarlo::attach_observer_to_lowest_score_conformation(
// 	core::conformation::ConformationObserverOP observer
// )
// {
// 	last_accepted_pose_->conformation().attach_observer( observer );
// }

// void
// DualMonteCarlo::attach_observer_to_last_accepted_pose(
// 	core::pose::PoseObserver * observer
// )
// {
// 	last_accepted_pose_->attach_observer( observer );
// }

// void
// DualMonteCarlo::attach_observer_to_lowest_score_pose( core::pose::PoseObserver * observer )
// {
// 	last_accepted_pose_->attach_observer( observer );
// }

core::scoring::ScoreFunction const &
DualMonteCarlo::score_function() const
{
	return *score_function_;
}

core::Real
DualMonteCarlo::last_accepted_score() const
{
	return last_accepted_pose_->energies().total_energy();
}


core::Real
DualMonteCarlo::lowest_score() const
{
	return lowest_score_pose_->energies().total_energy();
}

/// @brief const access to MC object
protocols::moves::MonteCarlo const &
DualMonteCarlo::MC() const
{ return MC_; }

/////////////////////////////////////////////////////////////////////////////
// void
// DualMonteCarlo::reset_counters()
// {
// 	trial_counter.clear();
// 	accept_counter.clear();
// 	energy_drop_counter.clear();
// }

// ///@detail return number of trials since last reset
// core::Size
// DualMonteCarlo::total_trials() const {
// 	Size ntrials( 0 );
// 	for ( std::map< std::string, int >::const_iterator
// 					it=trial_counter.begin(); it != trial_counter.end(); ++it ) {
// 		ntrials += it->second;
// 	}
// 	return ntrials;
// }
// /////////////////////////////////////////////////////////////////////////////
// void
// DualMonteCarlo::show_counters() const
// {
// 	for ( std::map< std::string, int >::const_iterator
// 					it=trial_counter.begin(); it != trial_counter.end(); ++it ) {
// 		std::string const & tag( it->first );
// 		int const ntrials( it->second );
// 		if  (accept_counter.count( tag )){
// 			int const accepts( accept_counter.find( tag )->second );
// 			core::Real const energy_drop( energy_drop_counter.find( tag )->second );
// 			TR << A( 16, tag ) <<
// 				" trials= " << I( 6, ntrials ) << "; " <<
// 				" accepts= " << F( 6, 4, core::Real( accepts )/ntrials ) << "; " <<
// 				" energy_drop/trial= " << F( 9, 5, core::Real( energy_drop ) / ntrials )<<
// 				std::endl;
// 		} else {
// 			TR << A( 16, tag ) <<
// 				" trials= " << I( 6, ntrials ) <<
// 				" NO ACCEPTS." << std::endl;
// 		} // else
// 	} // for
// }



} // namespace moves
} // namespace protocols
