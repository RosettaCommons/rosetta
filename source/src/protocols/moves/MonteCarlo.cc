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

static basic::Tracer TR("protocols.moves.MonteCarlo");
static numeric::random::RandomGenerator mc_RG(62452); // <- Magic number, do not change it!!!

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
	score_function_( src.score_function_ ? src.score_function_->clone() : static_cast< core::scoring::ScoreFunction * > (0) ),
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
	counter_( new TrialCounter ),
	update_boinc_( true ),
	total_score_of_last_considered_pose_( 0.0 ),
	last_accepted_score_( 0.0 ),
	lowest_score_( 0.0 ),
	heat_after_cycles_( 150 )
{
	last_accepted_pose_ = new Pose();
	lowest_score_pose_ = new Pose();
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
	counter_( new TrialCounter ),
	update_boinc_( true ),
	total_score_of_last_considered_pose_( 0.0 ),
	last_accepted_score_( 0.0 ),
	lowest_score_( 0.0 ),
	heat_after_cycles_( 150 )
{
	last_accepted_pose_ = new Pose();
	lowest_score_pose_ = new Pose();
	// score_function_ = new ScoreFunction(scorefxn);
	score_function_ = scorefxn.clone();
	last_check_ = 0;
	check_frequency_ = basic::options::option[ basic::options::OptionKeys::mc::convergence_check_frequency ]();
}


MonteCarlo::~MonteCarlo()
{
}


void MonteCarlo::clear_poses() {
	last_accepted_pose_ = new Pose();
	lowest_score_pose_ = new Pose();
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

void
MonteCarlo::set_counter(TrialCounterOP counter)
{
	counter_ = counter;
}

void
MonteCarlo::count_trial(std::string const & tag)
{
	counter_->count_trial(tag);
}

void
MonteCarlo::count_accepted(std::string const & tag)
{
	counter_->count_accepted(tag);
}

void
MonteCarlo::count_energy_drop(std::string const & tag, Real drop)
{
	counter_->count_energy_drop(tag, drop);
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

/////////////////////////////////////////////////////////////////////////////
TrialCounterCOP
MonteCarlo::counter() const {
	return counter_;
}

///@detail return number of trials since last reset
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
// behavior depends on setting of temperature
//
// return true for an accept, false otherwise
//
// 	mc_accepted
// 		3 = accepted:score beat low score and last_accepted score
// 		2 = accepted:score beat last_accepted score
// 		1 = thermally accepted: score worse than last_accepted score
// 		0 = not accepted

bool
MonteCarlo::boltzmann(
	Pose & pose,
	std::string const & move_type, // = unk
	core::Real const proposal_density_ratio, // = 1
	core::Real const inner_score_temperature_delta // = 0
)
{

// Work around a current bug in the pose observer classes..
#ifdef BOINC_GRAPHICS
	  if( update_boinc_ )
		boinc::Boinc::update_graphics_current( pose );
#endif


	// score the pose:
	Real const score( (*score_function_)( pose ) );
	total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.

	if ( false ) { // DEBUG
		pose::Pose copy_pose;
		copy_pose = pose;
		copy_pose.energies().clear();
		TR << "score copy" << std::endl;
		Real const copy_score = (*score_function_)( copy_pose );
		if ( std::abs( copy_score - score ) > 1E-6) {
			TR << "Score discrepancy.  score: " << score << " vs copy score: " << copy_score << std::endl;
			TR << "pose score: ";
			pose.energies().total_energies().show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;

			TR << "copy score: ";
			copy_pose.energies().total_energies().show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;

			TR << "Difference: ";
			core::scoring::EnergyMap emap = pose.energies().total_energies();
			emap -= copy_pose.energies().total_energies();
			emap.show_if_nonzero_weight( TR, score_function_->weights() );
			TR << std::endl;
		}
	}

	counter_->count_trial( move_type );


#ifdef BOINC_GRAPHICS
	if( update_boinc_ )
		boinc::Boinc::update_mc_trial_info( counter_->trial( move_type ), move_type );
#endif

	Real const boltz_factor = ( last_accepted_score() - score ) / temperature_ + inner_score_temperature_delta;
	Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) ) * proposal_density_ratio;
	if ( probability < 1 ) {
		if ( mc_RG.uniform() >= probability ) {
			mc_accepted_ = MCA_rejected; // rejected
			autotemp_reject();
			evaluate_convergence_checks( pose, true /*reject*/, false /* not final*/ );
			pose = ( *last_accepted_pose_ );

			//TR << "Rejected with score = " << score
			//	<< " ( prob = " << probability << " ) last_accepted scrore: " << last_accepted_score() << std::endl;
			// score_function_->show( T("protocols.moves.MonteCarlo.boltzmann"), pose );

			return false; // rejected


		}
		mc_accepted_ = MCA_accepted_thermally; // accepted thermally

		//TR  << "thermally accepted with score = " << score
		//	<< " ( prob = " << probability << " ) last_accepted scrore: " << last_accepted_score() << std::endl;
		// score_function_->show( T("protocols.moves.MonteCarlo.boltzmann"), pose );
	} else {

		mc_accepted_ = MCA_accepted_score_beat_last; // energy is lower than last_accepted
	}

	// these are useful but cost a little time to get
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( pose );

	counter_->count_accepted( move_type );
	counter_->count_energy_drop( move_type, score - last_accepted_score() );
	PROF_START( basic::MC_ACCEPT );
	*last_accepted_pose_ = pose;
	last_accepted_score_ = score;

#ifdef BOINC_GRAPHICS
	if( update_boinc_ )
		boinc::Boinc::update_graphics_last_accepted( pose, last_accepted_score() );
#endif

// 	last_accepted_pose_ = new Pose( *pose );

	autotemp_accept();

	if ( score < lowest_score() ) {
		*lowest_score_pose_ = pose;
		lowest_score_ = score;
		evaluate_convergence_checks( pose, false /*not reject*/, false /*not final*/ );

#ifdef BOINC_GRAPHICS
	if( update_boinc_ )
		boinc::Boinc::update_graphics_low_energy( pose, lowest_score() );
#endif

		//lowest_score_pose_ = new Pose ( pose );
		mc_accepted_ = MCA_accepted_score_beat_low; //3;
	}
	PROF_STOP( basic::MC_ACCEPT );
	return true; // accept!
}


bool
MonteCarlo::boltzmann(
	core::Real score_delta,
	std::string const & move_type, // = "unk"
	core::Real const proposal_density_ratio // = 1
)
{
	// score the pose:
	Real const score( last_accepted_score_ + score_delta );
	total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.

	counter_->count_trial( move_type );

	Real const boltz_factor = ( last_accepted_score() - score ) / temperature_;
	Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) ) * proposal_density_ratio;
	if ( probability < 1 ) {
		if ( mc_RG.uniform() >= probability ) {
			mc_accepted_ = MCA_rejected; // rejected
			autotemp_reject();

			//TR << "Rejected with score = " << score
			//	<< " ( prob = " << probability << " ) last_accepted scrore: " << last_accepted_score() << std::endl;
			// score_function_->show( T("protocols.moves.MonteCarlo.boltzmann"), pose );

			return false; // rejected
		}
		mc_accepted_ = MCA_accepted_thermally; // accepted thermally

		//TR  << "thermally accepted with score = " << score
		//	<< " ( prob = " << probability << " ) last_accepted scrore: " << last_accepted_score() << std::endl;
		// score_function_->show( T("protocols.moves.MonteCarlo.boltzmann"), pose );
	} else {

		mc_accepted_ = MCA_accepted_score_beat_last; // energy is lower than last_accepted
	}

	counter_->count_accepted( move_type );
	counter_->count_energy_drop( move_type, score - last_accepted_score() );
	PROF_START( basic::MC_ACCEPT );
	last_accepted_score_ = score;

	autotemp_accept();

	if ( score < lowest_score() ) {
		lowest_score_ = score;
		mc_accepted_ = MCA_accepted_score_beat_low; //3;
	}
	PROF_STOP( basic::MC_ACCEPT );
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


void
MonteCarlo::reset_last_accepted( Pose const & pose )
{
	PROF_START( basic::MC_ACCEPT );
	*last_accepted_pose_ = pose;
	PROF_STOP( basic::MC_ACCEPT );

	Real const score( (*score_function_)( *last_accepted_pose_ ) );
	/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );
	last_accepted_score_ = score;

	// if the last accepted pose has a lower energy that the lowest score pose, set that as well
	if (last_accepted_pose_->energies().total_energy() < lowest_score_pose_->energies().total_energy()) {
		PROF_START( basic::MC_ACCEPT );
		*lowest_score_pose_ = *last_accepted_pose_;
		PROF_STOP( basic::MC_ACCEPT );
	}
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


// for recovering from a checkpoint
void
MonteCarlo::set_last_accepted_pose( Pose const & pose )
{
	*last_accepted_pose_ = pose;
	last_accepted_score_ = last_accepted_pose_->energies().total_energy();
}

// for recovering from a checkpoint
void
MonteCarlo::set_lowest_score_pose( Pose const & pose )
{
	*lowest_score_pose_ = pose;
	lowest_score_ = lowest_score_pose_->energies().total_energy();
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
		temperature_ = std::min( temperature_ + heat_delta, max_temperature );
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
		temperature_ = quench_temp_;
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
