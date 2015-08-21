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
/// @author


#ifndef INCLUDED_protocols_moves_MonteCarlo_hh
#define INCLUDED_protocols_moves_MonteCarlo_hh


// type headers
#include <core/types.hh>

// unit headers
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/TrialCounter.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <protocols/moves/MonteCarloExceptionConverge.fwd.hh>
#include <utility/vector1.hh>
#include <string>

#ifdef WIN32
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#endif


// Forward declarations

namespace protocols {
namespace moves {

/// @brief This object is responsible for all of the major functions needed in
/// a Monte Carlo simulation. Its main purpose is to apply the Metropolis
/// Criterion on a pose, based on a ScoreFunction, temperature, and the
/// previously accepted pose. It stores the lowest-energy pose ecountered,
/// the last-accepted pose in the simulation, and various other statistics.
///
///
/// Output Methods:
///     MonteCarlo.show_counters()
///     MonteCarlo.show_scores()
///     MonteCarlo.show_state()
/// Common Methods:
///     MonteCarlo.last_accepted_score
///     MonteCarlo.last_accepted_pose
///     MonteCarlo.lowest_score
///     MonteCarlo.lowest_score_pose
///     MonteCarlo.score_function
///     MonteCarlo.set_temperature
///     MonteCarlo.temperature
class MonteCarlo : public utility::pointer::ReferenceCount {
public:
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::pose::PoseCOP PoseCOP;
	typedef core::Real Real;


public:

	/// @brief Copy constructor
	MonteCarlo( MonteCarlo const & );

	/// @brief Constructs a useable MonteCarlo object
	///
	/// mc = MonteCarlo( init_pose , scorefxn , temp )
	///
	/// Pose           init_pose   /manipulated during the simulation
	/// ScoreFunction  scorefxn    /evaluates pose scores
	/// Real (float)   temp        /used in the Metropolis Criterion
	MonteCarlo(
		Pose const & init_pose, // PoseCOP init_pose,
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature
	);


	/// @brief Constructor without Pose -- call reset(pose) before first use
	MonteCarlo(
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		Real const temperature
	);

	/// @brief Empty destructor in C++ file to reduce number of necessary includes
	virtual ~MonteCarlo();

	virtual
	MonteCarloOP clone() {
		return MonteCarloOP( new MonteCarlo( *this ) );
	}

	/// @brief Resets the ScoreFunction
	void
	reset_scorefxn(
		Pose const & init_pose,
		ScoreFunction const & scorefxn
	);

	/// @brief Change the weight on a score term in the object's scorefunction. Useful when we don't want to reset the whole scorefunction during an annealing step.
	void
	change_weight( core::scoring::ScoreType const & t, Real const & setting );

	/// @brief Sets the temperature value used in the Metropolis Criterion to  <temp>
	///
	/// example(s):
	///     mc.set_temperature( temp )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.temperature
	///     MonteCarlo.show_state
	virtual void
	set_temperature( Real const temp );


	/// @brief Returns the temperature value used in the Metropolis Criterion
	///
	/// example(s):
	///     mc.temperature()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.set_temperature
	///     MonteCarlo.show_state
	Real
	temperature() const {
		return temperature_;
	}

	/// @brief Sets autotemp to quench_temp
	/// example(s):
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.autotemp
	///     MonteCarlo.show_state
	void
	set_autotemp(
		bool const setting,
		core::Real const quench_temp
	);


	/// @brief Applies the Metropolis Criterion on pose based on
	/// the ScoreFunction, temperature, and the last accepted
	/// pose. This method evaluates the change in score, compares
	/// the trial pose to the last accepted pose, and updates the
	/// pose structure and simulation statistics appropriately
	///
	/// example(s):
	///     mc.boltzmann( pose )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	virtual bool
	boltzmann(
		Pose & pose,//PoseOP pose,
		std::string const & move_type = "unk",
		core::Real const proposal_density_ratio = 1,
		core::Real const inner_score_delta_over_temperature = 0
	);


	/// @brief Applies the Metropolis Criterion on the inputted
	/// pose based on the supplied score delta
	///
	/// example(s):
	///
	/// See also:
	///     MonteCarlo
	virtual bool
	boltzmann(
		core::Real score,
		std::string const & move_type = "unk",
		core::Real const proposal_density_ratio = 1,
		core::Real const inner_score_delta_over_temperature = 0,
		bool check_lowest_score = true
	);


	/// @brief Sets lowest score pose and last accepted pose to
	/// the score of  <pose>
	/// @note (does not reset counters)
	///
	/// example(s):
	///     mc.reset(pose)
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	///     MonteCarlo.lowest_scored_pose
	virtual
	void
	reset( Pose const & pose );

	/// @brief Sets the last accepted pose to the score of  <pose>
	/// @note (does not reset counters)
	void
	set_last_accepted_pose(
		Pose const & pose
	);

	/// @brief Sets the last accepted pose and last accepted score
	/// @note (does not reset counters)
	void
	set_last_accepted_pose(
		core::pose::Pose const& pose,
		core::Real score
	);

	void
	set_lowest_score_pose(
		core::pose::Pose const& pose
	);

	void
	set_lowest_score_pose(
		core::pose::Pose const& pose,
		core::Real score
	);

	/// @brief Returns the last accepted pose
	///
	/// example(s):
	///     mc.last_accepted_pose()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accept
	///     MonteCarlo.last_accepted_score
	Pose const &
	last_accepted_pose() const
	{
		return *last_accepted_pose_;
	}


	/// @brief Returns the lowest score pose encountered
	///
	/// example(s):
	///     mc.lowest_score_pose()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.lowest_score
	Pose const &
	lowest_score_pose() const
	{
		return *lowest_score_pose_;
	}


	/// @brief Compares score of <pose> to the lowest score found.
	/// If lower, sets the current lowest score pose and lowest score.
	/// Use internal pose energies if score_pose is false.
	/// Used to evaluate lowest score without boltzmann. Does not change pose structure.
	/// @note Does not update simulation statistics or last accepts by default.
	///
	/// example(s):
	///  mc.eval_lowest_score_pose( pose )
	/// See also:
	///  MonteCarlo
	///  MonteCarlo.lowest_score
	///  MonteCarlo.lowest_score_pose
	///  MonteCarlo.recover_low
	bool eval_lowest_score_pose(
		Pose & pose,
		bool score_pose = true,
		bool update_stats = false,
		std::string const & move_type = "unk"
	);


	/// @brief attach observer to last accepted conformation
	/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
	template< typename ConformationObserver >
	void
	attach_observer_to_last_accepted_conformation( ConformationObserver & obs );


	/// @brief attach observer to lowest score conformation
	/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
	template< typename ConformationObserver >
	void
	attach_observer_to_lowest_score_conformation( ConformationObserver & obs );


	/// @brief attach observer to last accepted pose
	/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
	template< typename PoseObserver >
	void
	attach_observer_to_last_accepted_pose( PoseObserver & obs );


	/// @brief attach observer to lowest score pose
	/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
	template< typename PoseObserver >
	void
	attach_observer_to_lowest_score_pose( PoseObserver & obs );


	/// @brief Sets the input  <pose>  and last accepted pose to
	/// the lowest score pose
	///
	/// example(s):
	///     mc.recover_low( pose )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accept
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	///     MonteCarlo.lowest_score_pose
	void
	recover_low( Pose & pose );


	/// @brief Sets the ScoreFunction to  <scorefxn> , re-scores
	/// last accepted pose and lowest score pose
	///
	/// example(s):
	///     mc.score_function( scorefxn )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.boltzmann
	///     MonteCarlo.set_temperature
	///     MonteCarlo.temperature
	///     ScoreFunction
	///     create_score_function
	virtual
	void
	score_function( ScoreFunction const & scorefxn ); // ScoreFunctionCOP scorefxn )


	/// @brief Returns the MonteCarlo ScoreFunction
	///
	/// example(s):
	///     mc.score_function()
	///     mc.score_function()( pose )
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.boltzmann
	///     MonteCarlo.set_temperature
	///     MonteCarlo.temperature
	///     ScoreFunction
	///     create_score_function
	ScoreFunction const & score_function() const;


	/// @brief Displays the last accepted score and the lowest score
	///
	/// example(s):
	///     mc.show_scores()
	/// Output as:
	///     protocols.moves.MonteCarlo: MonteCarlo:: last_accepted_score,lowest_score: X Y
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	///     MonteCarlo.show_counters
	///     MonteCarlo.show_state
	void show_scores() const;


	/// @brief Resets the mover counters
	///
	/// example(s):
	///     mc.reset_counters()
	/// See alse:
	///     MonteCarlo
	///     MonteCarlo.show_counters
	void reset_counters();


	/// @brief Displays the entire MonteCarlo state
	/// temperature, scores, annealing settings,
	/// move statistics, move counters (show_counters)
	///
	/// example(s):
	///     mc.show_state()
	/// Output as:
	///     protocols.moves.MonteCarlo: MC: t l1 l2 las lws la au qu mca
	///         t= temperature
	///         l1= (*score_function_)(*last_accepted_pose_)
	///         l2= (*score_function_)(*lowest_score_pose_)
	///         las= last accepted score
	///         lws= lowest score
	///         la= last_accept_
	///         au= autotemp_
	///         qu= quench_temp_
	///         mca= mc_accepted_
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.show_counters
	///     MonteCarlo.show_scores
	///     MonteCarlo.last_accepted_score
	///     MonteCarlo.lowest_score
	///     MonteCarlo.temperature
	void show_state() const;


	/// @brief Displays the number of trials performed, fraction
	/// of trial moves accepted, and the average energy drop per
	/// accepted trial by mover types applied (unknown movers or
	/// perturbations are listed as "unktrials")
	///
	/// example(s):
	///     mc.show_counters()
	/// Output as:
	///     protocols.moves.MonteCarlo:            unk trials=     X;  accepts=     Y;  energy_drop/trial=     Z
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.show_scores
	///     MonteCarlo.show_state
	void show_counters() const;

	/// @brief Set the counter to be used for this simulation.
	// Undefined, commenting out to fix PyRosetta build  void set_counter(TrialCounterOP counter);

	/// @brief Return statistics for each type of move that has been attempted
	/// during this simulation.
	///
	/// example(s):
	///     mc.counters()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.show_counters
	///     MonteCarlo.reset_counters
	// Undefined, commenting out to fix PyRosetta build  TrialCounterCOP counter() const;

	/// @brief Returns the total number of trials since the last reset
	/// @note: MonteCarlo.boltzmann(pose) updates the number of trials
	///
	/// example(s):
	///     mc.total_trials()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accept
	///     MonteCarlo.show_counters
	///     MonteCarlo.show_state
	Size total_trials() const;

	/// @brief Manually record that a move was attempted.
	/// @details This is intended for moves that are made outside to context of
	/// MonteCarlo itself.  Temperature moves are a good example.
	// Undefined, commenting out to fix PyRosetta build  void count_trial(std::string const & tag);

	/// @brief Manually record that a move was accepted.
	/// @details This is intended for moves that are made outside to context of
	/// MonteCarlo itself.  Temperature moves are a good example.
	// Undefined, commenting out to fix PyRosetta build  void count_accepted(std::string const & tag);

	/// @brief Manually record that a move produced a change in energy.
	/// @details This is intended for moves that are made outside to context of
	/// MonteCarlo itself.  Temperature moves are a good example.
	// Undefined, commenting out to fix PyRosetta build void count_energy_drop(std::string const & tag, core::Real drop);

	/// @brief Returns the score value of the last accepted pose
	///
	/// example(s):
	///     mc.last_accepted_score()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accept
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.show_counters
	///     MonteCarlo.show_scores
	///     MonteCarlo.show_state
	Real last_accepted_score() const;


	/// @brief Returns the score value of the lowest score pose encountered
	///
	/// example(s):
	///     mc.lowest_score()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.lowest_score_pose
	///     MonteCarlo.show_counters
	///     MonteCarlo.show_scores
	///     MonteCarlo.show_state
	Real lowest_score() const;


	/// @brief Returns mc_accepted, informative of the last move applied
	///
	/// Note: Returns true for an accept, false otherwise
	///     3 = accepted:score beat low score and last_accepted score
	///     2 = accepted:score beat last_accepted score
	///     1 = thermally accepted: score worse than last_accepted score
	///     0 = not accepted
	/// example(s):
	///     mc.mc_accepted()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.show_state
	MCA mc_accepted() const;

	/// @brief Returns mc_accepted as a human-readable string, informative of the last move applied
	std::string mc_accepted_string() const;

	/// @brief Removes last accepted pose and lowest score pose
	///
	/// example(s):
	///     mc.clear_poses()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.lowest_score_pose
	///     MonteCarlo.recover_low
	///     MonteCarlo.reset
	///     MonteCarlo.set_last_accepted_pose
	///     MonteCarlo.set_lowest_score_pose
	void clear_poses();// remove last_accepted and lowest_score


	/// no brief for now
	void set_update_boinc( bool setting ){ update_boinc_ = setting; }
	bool get_update_boinc() const { return update_boinc_; }


	Real total_score_of_last_considered_pose() const { return total_score_of_last_considered_pose_; }

	/// @brief Returns the number of trials since last acceptance
	///
	/// example(s):
	///     mc.last_accept()
	/// See also:
	///     MonteCarlo
	///     MonteCarlo.show_counters
	///     MonteCarlo.last_accepted_pose
	///     MonteCarlo.last_accepted_score
	core::Size last_accept() const {
		return last_accept_;
	}


	/// no brief for now
	core::Size heat_after_cycles() const {
		return heat_after_cycles_;
	}


	/// no brief for now
	void set_heat_after_cycles( core::Size setting ) {
		heat_after_cycles_ = setting;
	}


	/// no brief for now
	void push_back( moves::MonteCarloExceptionConvergeOP );
	/////////////////////////////////////////////////////////////////////////////
	// private methods
protected:

	/// @brief for managing the temperature, if we need to do so
	void
	autotemp_reject();

	void
	autotemp_accept();

	void
	evaluate_convergence_checks( core::pose::Pose const& pose, bool reject, bool final );

public:
	Size
	check_frequency() const {
		return check_frequency_;
	}

private:
	/// unimplemented -- do not use
	MonteCarlo const & operator = ( MonteCarlo const & ); // assignment operator -- do not use.


protected:
	void set_last_accepted_score(
		core::Real score
	) {
		last_accepted_score_ = score;
	};

	void
	set_mc_accepted( MCA value ) {
		mc_accepted_ = value;
	}
	/////////////////////////////////////////////////////////////////////////////
	// data
private:

	/// @brief Latest accepted pose
	PoseOP last_accepted_pose_;

	/// @brief Lowest score pose encountered
	PoseOP lowest_score_pose_;

	/// @brief Acceptance criterion temperature
	core::Real temperature_;

	/// @brief Internal scoring function
	ScoreFunctionOP score_function_;

	/// @brief For abinitio-style increasing the temperature after a large number of rejects:
	bool autotemp_;
	core::Real quench_temp_;
	int last_accept_;

	/// @brief Result of the last call to boltzmann
	MCA mc_accepted_;

	/// @brief diagnostics
	TrialCounterOP counter_;

	bool update_boinc_;

	Real total_score_of_last_considered_pose_;
	Real last_accepted_score_;
	Real lowest_score_;

	Size heat_after_cycles_;

	utility::vector1< moves::MonteCarloExceptionConvergeOP > convergence_checks_;
	Size last_check_;
	Size check_frequency_;
};

// for Python bindings
std::ostream & operator << ( std::ostream & os, MonteCarlo const & mc);

} // moves
} // rosetta

#endif
