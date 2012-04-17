// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TrialMover
/// @brief performs a move and accepts it according to Monte Carlo accept/reject criterion.
/// @author Monica Berrondo

#include <protocols/moves/TrialMover.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.hh>

// Utility headers

// C++ headers
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

basic::Tracer tr("protocols.TrialMover");
	
using namespace core;

MonteCarloUtil::MonteCarloUtil() : mc_(0)
{
	
}
	
MonteCarloUtil::MonteCarloUtil( protocols::moves::MonteCarloOP mc) : mc_(mc)
{
	
}
	
	MonteCarloUtil::~MonteCarloUtil() {}
	
void MonteCarloUtil::apply(Pose & pose)
{
	if(mode_ == "reset")
	{
		mc_->reset(pose);
	}else if(mode_ == "recover_low")
	{
		mc_->recover_low(pose);
	}else
	{
		utility_exit_with_message("MonteCarloUtil mode must be 'reset' or 'recover_low', this should have been caught earlier, dying");
	}
	
}

void MonteCarloUtil::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & /* pose */)
{
	if(!tag->hasOption("mode"))
	{
		utility_exit_with_message("you must specify option mode in MonteCarloUtil");
	}
	if(!tag->hasOption("montecarlo"))
	{
		utility_exit_with_message("you must specify the option montecarlo in MonteCarloUtil");
	}
	
	mode_ = tag->getOption<std::string>("mode");
	std::string const mc_name(tag->getOption<std::string>("montecarlo"));
	
	if(mode_ != "reset" && mode_ != "recover_low")
	{
		utility_exit_with_message("the option mode must be set to either reset or recover_low in MonteCarloUtil");
	}
	
	mc_ = *data.get<protocols::moves::MonteCarlo*>( "montecarlos",mc_name);
}
	
protocols::moves::MoverOP MonteCarloUtil::clone() const
{
	return new MonteCarloUtil(mc_);
}



std::string MonteCarloUtil::get_name() const
{
	return "MonteCarloUtil";
}

	
TrialMover::TrialMover() :
	start_weight_( 0.0 ),
	original_weight( 0.0 ),
	ramp_weight( 0.0 ),
	delta( 0.0 ),
	stats_type_( all_stats )
{}

// constructor with arguments
TrialMover::TrialMover( MoverOP mover_in, MonteCarloOP mc_in ) :
	start_weight_( 0.0 ),
	original_weight( 0.0 ),
	ramp_weight( 0.0 ),
	delta( 0.0 ),
	stats_type_( all_stats )
{
	mover_ = mover_in;
	mc_ = mc_in;
}

TrialMover::~TrialMover() {}

// set the weights for the score_type for ramping
void TrialMover::initialize_weights(
	Real const start_weight,
	Real const end_weight,
	core::scoring::ScoreType score_type,
	int const ramp_cycles
) {
	original_weight = mc_->score_function().get_weight( score_type );
	delta = ( end_weight - start_weight ) / ramp_cycles;
	ramp_weight = start_weight;
	start_weight_ = start_weight;
}

void TrialMover::set_mc( MonteCarloOP mc_in ) {
	mc_ = mc_in;
}


/// @begin TrialMover::apply()
/// @brief:
/// 	the apply function for a trial
///	@detailed:
///		the trial object is created with an mc object
///		the mover is applied before doing an mc.boltzmann
///
///	@author: Monica Berrondo
void TrialMover::apply( pose::Pose & pose )
{	
	using scoring::total_score;

	/// get the initial scores
	if ( keep_stats_type() == all_stats ) {
		stats_.add_score( mc_->last_accepted_score() ); ///< initial_last_accepted_score
		stats_.add_score( pose.energies().total_energy() ); ///< initial_pose_score
	}

	/// make the move
	mover_->apply( pose );

	// if ( keep_stats_type() == all_stats ) { //// score and get residue energies
	// Stupid and wasteful.  The structure will be scored inside mc_->boltzman.  mc_->score_function()( pose );
	// Unneccessary since revision 23846 --- mc_->score_function().accumulate_residue_total_energies( pose );
	// WAIT FOR IT. stats_.add_score( pose.energies().total_energy() ); ///< score_after_move
	// }

	/// test if MC accepts or rejects it
	bool accepted_move = mc_->boltzmann( pose, mover_->type() );

	if ( keep_stats_type() == all_stats ) {
		stats_.add_score( mc_->total_score_of_last_considered_pose() );
	}

	if ( stats_type_ <= accept_reject ) {
		stats_.accepted( accepted_move );
	}
	if ( keep_stats_type() > no_stats ) {
		stats_.print( mc_, mover_->type() );
	}
}

std::string
TrialMover::get_name() const {
	return "TrialMover";
}


Real TrialMover::acceptance_rate() const
{
	tr << "Acceptance rate: " << stats_.acceptance_rate() << std::endl;
	return stats_.acceptance_rate();
}

int TrialMover::num_accepts() const
{
	return stats_.num_accepted();
}

// sets the input pose also for the contained mover (barak)
void TrialMover::set_input_pose( PoseCOP pose )
{
	this->Mover::set_input_pose( pose );
	if(mover_)
		mover_->set_input_pose( pose );
}

// sets the native pose also for the contained mover (barak)
void TrialMover::set_native_pose( PoseCOP pose )
{
	this->Mover::set_native_pose( pose );
	if(mover_)
		mover_->set_native_pose( pose );
}

void TrialMover::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const &,
	Movers_map const & movers,
	Pose const &
)
{
	// 1. MonteCarlo object's name
	std::string const mc_name( tag->getOption< std::string > ( "montecarlo", "" ));
	if ( mc_name == "" ) {
		utility_exit_with_message( "TrialMover requires the 'montecarlo' option which was not provided" );
	}

	// 2. Mover
	std::string const movername( tag->getOption< std::string > ( "mover", "" ));
	if ( movername == "" ) {
		utility_exit_with_message( "TrialMover requires the 'mover' option which was not provided" );	
	}
	Movers_map::const_iterator  find_mover ( movers.find( movername ));
	if( find_mover == movers.end() && movername != "" ) {
		utility_exit_with_message( "TrialMover was not able to find the mover named '" + movername + "' in the Movers_map" );
	}
	mover_ = find_mover->second;
	mc_ = *data.get< protocols::moves::MonteCarlo * >( "montecarlos", mc_name );

	// 3. stats_type.
	std::string const statstype( tag->getOption< std::string > ( "keep_stats", "no_stats" ));
	if ( statstype != "no_stats" && statstype != "accept_reject" && statstype != "all_stats" ) {
		utility_exit_with_message( "TrialMover keep_stats may only be given the values:\n'no_stats', 'accept_reject', and 'all_stats'.\nRead value '" + statstype + "'" );
	}
	if ( statstype == "no_stats" ) {
		stats_type_ = no_stats;
	} else if ( statstype == "accept_reject" ) {
		stats_type_ = accept_reject;
	} else {
		stats_type_ = all_stats;
	}

}


}  // namespace moves
}  // namespace protocols
