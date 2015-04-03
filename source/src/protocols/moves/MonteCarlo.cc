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
/// @author Phil Bradley

// Unit Headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialCounter.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> //Pretty output.

// Numeric Headers
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>

#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <utility/vector1.hh>


// for rosetta++ like boinc graphics
#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.moves.MonteCarlo" );

namespace protocols {
namespace moves {

using namespace core;
using namespace ObjexxFCL::format;

/// @details The copy constructor does not copy the OPs, but rather creates new objects using the copy
/// constructors or the clone() methods of the objects being pointed at.  This is important, since otherwise,
/// a copy of a Monte Carlo object could corrupt the state held in the original MonteCarlo object.
MonteCarlo::MonteCarlo( MonteCarlo const & src ) :
	utility::pointer::ReferenceCount(),
	last_accepted_pose_( src.last_accepted_pose_ ? new core::pose::Pose( * src.last_accepted_pose_ ) : 0 ),
	lowest_score_pose_( src.lowest_score_pose_ ? new core::pose::Pose( * src.lowest_score_pose_ ) : 0 ),
	temperature_( src.temperature_ ),
	score_function_( src.score_function_ ? src.score_function_->clone() : core::scoring::ScoreFunctionOP(0) ),
	autotemp_( src.autotemp_ ),
	quench_temp_( src.quench_temp_ ),
	last_accept_( src.last_accept_ ),
	mc_accepted_( src.mc_accepted_ ),
	counter_( src.counter_ ),
	update_boinc_( src.update_boinc_ ),
	total_score_of_last_considered_pose_( src.total_score_of_last_considered_pose_ ),
	last_accepted_score_( src.last_accepted_score_ ),
	lowest_score_( src.lowest_score_ ),
	heat_after_cycles_( src.heat_after_cycles_ ),
	convergence_checks_( src.convergence_checks_ ),
	last_check_( src.last_check_ ),
	check_frequency_( src.check_frequency_ )
{
}


// constructor for monte_carlo object
MonteCarlo::MonteCarlo(
	Pose const & init_pose, // PoseCOP init_pose,
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature
):
	temperature_( temperature ),
	autotemp_( false ),
	quench_temp_( 0.0 ),
	last_accept_( 0 ),
	mc_accepted_( MCA_accepted_score_beat_last ), // init_pose beats the absence of a pose
	counter_( TrialCounterOP( new TrialCounter ) ),
	update_boinc_( true ),
	total_score_of_last_considered_pose_( 0.0 ),
	last_accepted_score_( 0.0 ),
	lowest_score_( 0.0 ),
	heat_after_cycles_( 150 )
{
	last_accepted_pose_ = PoseOP( new Pose() );
	lowest_score_pose_ = PoseOP( new Pose() );
	// score_function_ = new ScoreFunction(scorefxn);
	score_function_ = scorefxn.clone();
	reset( init_pose );

	last_check_ = 0;
	check_frequency_ = basic::options::option[ basic::options::OptionKeys::mc::convergence_check_frequency ]();
}

MonteCarlo::MonteCarlo(
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature
):
	temperature_( temperature ),
	autotemp_( false ),
	quench_temp_( 0.0 ),
	last_accept_( 0 ),
	mc_accepted_( MCA_accepted_score_beat_last ), // an empty pose beats the absence of a pose
	counter_( TrialCounterOP( new TrialCounter ) ),
	update_boinc_( true ),
	total_score_of_last_considered_pose_( 0.0 ),
	last_accepted_score_( 0.0 ),
	lowest_score_( 0.0 ),
	heat_after_cycles_( 150 )
{
	last_accepted_pose_ = PoseOP( new Pose() );
	lowest_score_pose_ = PoseOP( new Pose() );
	// score_function_ = new ScoreFunction(scorefxn);
	score_function_ = scorefxn.clone();
	last_check_ = 0;
	check_frequency_ = basic::options::option[ basic::options::OptionKeys::mc::convergence_check_frequency ]();
}


MonteCarlo::~MonteCarlo()
{
}


void MonteCarlo::clear_poses() {
	last_accepted_pose_ = PoseOP( new Pose() );
	lowest_score_pose_ = PoseOP( new Pose() );
}

void
MonteCarlo::reset_scorefxn(
	Pose const & init_pose,
	ScoreFunction const & scorefxn
){
	// score_function_ = new ScoreFunction(scorefxn);
	score_function_ = scorefxn.clone();
	reset(init_pose);
}

void
MonteCarlo::set_temperature( Real const temp )
{
	temperature_ = temp;
}

/// return the simulation state to the lowest energy structure we've seen
void
MonteCarlo::recover_low( Pose & pose )
{
	( pose ) = ( *lowest_score_pose_ );
	*last_accepted_pose_ = *lowest_score_pose_ ;
	//last_accepted_pose_ = new Pose( *lowest_score_pose_ );
	last_accepted_score_ = last_accepted_pose_->energies().total_energy();
}


void
MonteCarlo::show_state() const
{
	TR << "MC: " << temperature_
	    << "  " << (*score_function_)(*last_accepted_pose_)
	    << "  " << (*score_function_)(*lowest_score_pose_)
			<< "  " << last_accepted_score()
			<< "  " << lowest_score()
	    << "  " << last_accept_
	    << "  " << autotemp_
	    << "  " << quench_temp_
		  << "  " << to_string( mc_accepted_ )
	    << std::endl;
	show_counters();
}

/////////////////////////////////////////////////////////////////////////////

void
MonteCarlo::show_scores() const
{
	TR << "MonteCarlo:: last_accepted_score,lowest_score: " <<
		last_accepted_score() << ' ' << lowest_score() << std::endl;
}
/////////////////////////////////////////////////////////////////////////////
void
MonteCarlo::reset_counters()
{
	counter_->reset();
}

/// @detail return number of trials since last reset
core::Size
MonteCarlo::total_trials() const {
	return counter_->total_trials();
}
/////////////////////////////////////////////////////////////////////////////
void
MonteCarlo::show_counters() const {
	counter_->show();
}

void
MonteCarlo::change_weight( core::scoring::ScoreType const & t, Real const & setting ) {
	score_function_->set_weight( t, setting );
}

/// set the scorefxn,  re-scores last-accepted and lowest-score pose
void
MonteCarlo::score_function( ScoreFunction const & scorefxn )
{
	using namespace scoring;

	// *score_function_ = scorefxn;
	// score_function_ = new ScoreFunction(scorefxn);
	score_function_ = scorefxn.clone();
	//TR << "new score_function_ within mc!" << std::endl;

	lowest_score_ = (*score_function_)( *lowest_score_pose_ );
	//TR << "lowest_score: " << lowest_score_ << " total: " << lowest_score_pose_->energies().total_energy() << std::endl;
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *lowest_score_pose_ );

	if ( false ) { // DEBUG
		pose::Pose copy_pose;
		copy_pose = *lowest_score_pose_;
		copy_pose.energies().clear();
		TR << "score copy" << std::endl;
		Real const copy_score = (*score_function_)( copy_pose );
		if ( std::abs( copy_score - lowest_score_ ) > 1E-6) {
			TR << "Score discrepancy.  lowest_score: " << lowest_score_ << " vs copy score: " << copy_score << std::endl;
			TR << "pose score: ";
			lowest_score_pose_->energies().total_energies().show_if_nonzero_weight( TR , score_function_->weights() );
			TR << std::endl;

			TR << "copy score: ";
			copy_pose.energies().total_energies().show_if_nonzero_weight( TR , score_function_->weights() );
			TR << std::endl;

			TR << "Difference: ";
			core::scoring::EnergyMap emap = lowest_score_pose_->energies().total_energies();
			emap -= copy_pose.energies().total_energies();
			emap.show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;
		}
	}

	//T("protocols.moves.MonteCarlo.score_function") << "lowest_score";
	//score_function_->show( T("protocols.moves.MonteCarlo.score_function"), *lowest_score_pose_ );
	//T("protocols.moves.MonteCarlo.score_function");

	last_accepted_score_ = (*score_function_)( *last_accepted_pose_ );
	//TR << "last_accepted_score: " << last_accepted_score << " total: " << last_accepted_pose_->energies().total_energy() << std::endl;
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );


	if ( false ) { // DEBUG
		pose::Pose copy_pose;
		copy_pose = *last_accepted_pose_;
		copy_pose.energies().clear();
		TR << "score copy" << std::endl;
		Real const copy_score = (*score_function_)( copy_pose );
		if ( std::abs( copy_score - last_accepted_score_ ) > 1E-6) {
			TR << "Score discrepancy.  last_accepted_score: " << last_accepted_score_ << " vs copy score: " << copy_score << std::endl;
			TR << "pose score: ";
			last_accepted_pose_->energies().total_energies().show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;

			TR << "copy score: ";
			copy_pose.energies().total_energies().show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;

			TR << "Difference: ";
			core::scoring::EnergyMap emap = last_accepted_pose_->energies().total_energies();
			emap -= copy_pose.energies().total_energies();
			emap.show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;
		}
	}


	//T("protocols.moves.MonteCarlo.score_function") << "last_accepted";
	//score_function_->show( T("protocols.moves.MonteCarlo.score_function"), *last_accepted_pose_ );
	//T("protocols.moves.MonteCarlo.score_function") << " finished!";

	///SML 6/25/08 - if last accepted is now better than lowest, change it!
	///explicit this-> necessary to resolve function name against identically named local variables
	if ( this->last_accepted_score() < this->lowest_score() ) {
		PROF_START( basic::MC_ACCEPT );
		*lowest_score_pose_ = *last_accepted_pose_;
		lowest_score_ = last_accepted_score_;
		PROF_STOP( basic::MC_ACCEPT );
	}

}


/////////////////////////////////////////////////////////////////////////////
//////////////////////////////
// @details behavior depends on setting of temperature
//
// return true for an accept, false otherwise
//
// 	mc_accepted
// 		3 = accepted:score beat low score and last_accepted score
// 		2 = accepted:score beat last_accepted score
// 		1 = thermally accepted: score worse than last_accepted score
// 		0 = not accepted
//
// Optional inputs:
//
//  proposal_density_ratio = ratio of backward proposal probability to
//                            forward proposal probability,
//                            to maintain detailed balance.
//
//  inner_score_delta_over_temperature
//                          = change in energy in any separate inner loop that
//                            used Boltzmann criterion with different energy
//                            function. Don't penalize with that energy
//                            difference again. See, e.g.,
//                             Hetenyi et al., JCP 117: 8203-8207.
//                            ( Note: functionally redundant with
//                               proposal_density_ratio)

bool
MonteCarlo::boltzmann(
	Pose & pose,
	std::string const & move_type, // = unk
	core::Real const proposal_density_ratio, // = 1
	core::Real const inner_score_delta_over_temperature // = 0
)
{

// Work around a current bug in the pose observer classes..
#ifdef BOINC_GRAPHICS
	  if( update_boinc_ )
		boinc::Boinc::update_graphics_current( pose );
#endif

	// score the pose:
	Real const score( (*score_function_)( pose ) );
	//now delegate decision making...
	bool const accept( boltzmann( score, move_type, proposal_density_ratio, inner_score_delta_over_temperature ) );

	//rejected ?
	if ( !accept ) {
		evaluate_convergence_checks( pose, true /*reject*/, false /* not final*/ );
		pose = ( *last_accepted_pose_ );
		return false; // rejected
	}

	//accepted !
	PROF_START( basic::MC_ACCEPT );
	*last_accepted_pose_ = pose;

	// print out the scores for each decoy to cmd out, if you pass a flag
	// nice for testing
	if ( basic::options::option[ basic::options::OptionKeys::mc::log_scores_in_MC ].user() && basic::options::option[ basic::options::OptionKeys::mc::log_scores_in_MC ]() == true ){
		score_function_->show( pose );
	}

#ifdef BOINC_GRAPHICS
	if( update_boinc_ )
		boinc::Boinc::update_graphics_last_accepted( pose, last_accepted_score() );
#endif

	if ( mc_accepted_ == MCA_accepted_score_beat_low ) {
		*lowest_score_pose_ = pose;
		evaluate_convergence_checks( pose, false /*not reject*/, false /*not final*/ );

#ifdef BOINC_GRAPHICS
		if( update_boinc_ )
			boinc::Boinc::update_graphics_low_energy( pose, lowest_score() );
#endif

	} //MCA_accepted_score_beat_low

	PROF_STOP( basic::MC_ACCEPT );
	return true; // accept!
}

//////////////////////////////////////////////////////////
// See notes above on boltzmann()
//////////////////////////////////////////////////////////
bool
MonteCarlo::boltzmann(
	core::Real score,
	std::string const & move_type, // = "unk"
	core::Real const proposal_density_ratio, // = 1
	core::Real const inner_score_delta_over_temperature, // = 0
	bool check_lowest_score //=true
)
{
	// figure out the actual score
	total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.
	counter_->count_trial( move_type );

#ifdef BOINC_GRAPHICS
	if( update_boinc_ )
		boinc::Boinc::update_mc_trial_info( counter_->trial( move_type ), move_type );
#endif

	Real const score_delta( score - last_accepted_score_ );
	Real const boltz_factor =  ( -score_delta / temperature_ ) + inner_score_delta_over_temperature;
	Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) ) * proposal_density_ratio;
	if ( probability < 1 ) {
		if ( numeric::random::rg().uniform() >= probability ) {
			mc_accepted_ = MCA_rejected; // rejected
			autotemp_reject();
			return false; // rejected
		}
		mc_accepted_ = MCA_accepted_thermally; // accepted thermally
	} else {
		mc_accepted_ = MCA_accepted_score_beat_last; // energy is lower than last_accepted
	}

	counter_->count_accepted( move_type );
	counter_->count_energy_drop( move_type, score_delta );
	last_accepted_score_ = score;

	autotemp_accept();

	if ( check_lowest_score && score < lowest_score() ) {
		lowest_score_ = score;
		mc_accepted_ = MCA_accepted_score_beat_low; //3;
	}
	return true; // accept!
}


void
MonteCarlo::reset( Pose const & pose )
{
	PROF_START( basic::MC_ACCEPT );
	*last_accepted_pose_ = pose;
	PROF_STOP( basic::MC_ACCEPT );

	Real const score( (*score_function_)( *last_accepted_pose_ ) );
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );

	PROF_START( basic::MC_ACCEPT );
	*lowest_score_pose_ = *last_accepted_pose_;
	PROF_STOP( basic::MC_ACCEPT );

	last_accepted_score_ = score;
	lowest_score_ = score;
}

/////////////////////////////////////////////////////////////////////////////
void
MonteCarlo::set_autotemp(
	bool const setting,
	Real const quench_temp
)
{
	autotemp_ = setting;
	quench_temp_ = quench_temp;
	last_accept_ = 0;
}

void
MonteCarlo::set_last_accepted_pose( Pose const & pose )
{
	*last_accepted_pose_ = pose;
	last_accepted_score_ = last_accepted_pose_->energies().total_energy();
}

void MonteCarlo::set_last_accepted_pose( core::pose::Pose const& pose, core::Real score ) {
	*last_accepted_pose_ = pose;
	last_accepted_score_ = score;
}

void MonteCarlo::set_lowest_score_pose( core::pose::Pose const& pose ) {
	*lowest_score_pose_ = pose;
	lowest_score_ = pose.energies().total_energy();
}

void MonteCarlo::set_lowest_score_pose( core::pose::Pose const& pose, core::Real score ) {
	*lowest_score_pose_ = pose;
	lowest_score_ = score;
}

bool
MonteCarlo::eval_lowest_score_pose(
	Pose & pose,
	bool score_pose, // true
	bool update_stats, //false
	std::string const & move_type //unk
)
{
	//Get or calculate energy
	Real score;
	if (score_pose){
		score = (*score_function_)(pose);
	}
	else{
		score =  pose.energies().total_energy();
	}

	//Evaluate
	total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.
	if ( score < lowest_score() ) {
		*lowest_score_pose_ = pose;
		lowest_score_ = score;
		if (update_stats){
			counter_->count_accepted( move_type );
			counter_->count_energy_drop( move_type, score - last_accepted_score() );
			last_accepted_score_ = score;
			mc_accepted_ = MCA_accepted_score_beat_low;
			*last_accepted_pose_ = pose;
			evaluate_convergence_checks( pose, false /*not reject*/, false /*not final*/ );
		}

		return true;
	}
	else{
		if (update_stats){
			mc_accepted_ = MCA_rejected; // rejected
		}
		return false;
	}
}

core::scoring::ScoreFunction const &
MonteCarlo::score_function() const
{
	return *score_function_;
}

Real
MonteCarlo::last_accepted_score() const
{
	//if (last_accepted_pose_->energies().total_energy() != last_accepted_score_) {
	//	TR << "Last: " << last_accepted_pose_->energies().total_energy() << " != " << last_accepted_score_ << " diff " << last_accepted_pose_->energies().total_energy() - last_accepted_score_ << std::endl;
	//}
	//return last_accepted_pose_->energies().total_energy();
	return last_accepted_score_;
}


Real
MonteCarlo::lowest_score() const
{
	//if (lowest_score_pose_->energies().total_energy() != lowest_score_) {
	//	TR << "Low: " << lowest_score_pose_->energies().total_energy() << " != " << lowest_score_ << " diff " << lowest_score_pose_->energies().total_energy() - lowest_score_ << std::endl;
	//}
	//return lowest_score_pose_->energies().total_energy();
	return lowest_score_;
}


MCA
MonteCarlo::mc_accepted() const
{
	return mc_accepted_;
}

std::string
MonteCarlo::mc_accepted_string() const
{
	return to_string( mc_accepted_ );
}


/////////////////////////////////////////////////////////////////////////////
// replicate logic from monte_carlo.cc
//
// should probably make these parameters ( 150, 1.0 )
void
MonteCarlo::autotemp_reject()
{
	//	int const heat_after_cycles( 150 );
	Real const heat_delta( quench_temp_ * 0.5 );
	Real const max_temperature( quench_temp_ * 10.0 );

	if ( !autotemp_ ) return;
	if ( last_accept_ >= (int) heat_after_cycles_ ) {
		//if ( temperature_ > max_temperature * 0.25 )
		TR << "autotemp_reject -- heat: " << last_accept_<< ' ' << temperature_  << std::endl;
		last_accept_ = -1;
		set_temperature( std::min( temperature_ + heat_delta, max_temperature ) );
	}
	++last_accept_;
}

/////////////////////////////////////////////////////////////////////////////
// replicate logic from monte_carlo.cc
void
MonteCarlo::autotemp_accept()
{
	if ( !autotemp_ ) return;
	if ( temperature_ != quench_temp_ ) {
		set_temperature( quench_temp_ );
		TR << "autotemp_accept: reset temperature_ = " << temperature_ << std::endl;
	}

	last_accept_ = 0;
}

void
MonteCarlo::evaluate_convergence_checks( core::pose::Pose const& pose, bool reject, bool /*final*/ ) {
	if ( !reject || numeric::mod( last_check_++, check_frequency_ ) == 0 ) {
		for ( utility::vector1< moves::MonteCarloExceptionConvergeOP >::iterator it = convergence_checks_.begin(); it != convergence_checks_.end(); ++it ) {
			(**it)( pose, *this, reject );
		}
	}
}

void
MonteCarlo::push_back( moves::MonteCarloExceptionConvergeOP check ) {
	convergence_checks_.push_back( check );
}

// for Python bindings
std::ostream & operator << ( std::ostream & os, MonteCarlo const & mc)
{
	os << "protocols.moves.MonteCarlo:\n";
	os << "\"Temperature\" (kT): " << mc.temperature() << "\n";
	os << "Total Trials: " << mc.total_trials() << "\n";
	os << "Lowest Score: " << mc.lowest_score() << "\n";
	os << "Last Accepted Score: " << mc.last_accepted_score() << "\n";
	os << "last_accept = " << mc.lowest_score() << std::endl;
	return os;
}

} // namespace moves
} // namespace core
