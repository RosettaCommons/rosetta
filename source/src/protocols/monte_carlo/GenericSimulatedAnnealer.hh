// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file src/protocols/monte_carlo/GenericSimulatedAnnealer.hh
/// @brief perform a given mover and sample structures by MonteCarlo with simulated annealing of acceptance temperatures
/// @details The "score" evaluation of pose during MC after applying mover is done by
/// ither FilterOP that can do report_sm() or ScoreFunctionOP you gave.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_GenericSimulatedAnnealer_hh
#define INCLUDED_protocols_simple_moves_GenericSimulatedAnnealer_hh

// Unit header
#include <protocols/monte_carlo/GenericSimulatedAnnealer.fwd.hh>
#include <basic/datacache/DataMapObj.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/monte_carlo/GenericMonteCarloMover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// External headers

namespace protocols {
namespace monte_carlo {

/// @brief Describes result of a trial
/// @details REJECTED : Move was rejected
///          ACCEPTED : Move was accepted
///          FAILED   : Some failure occurred during the move
///          FINISHED : The stopping condition was met
enum TrialResult {
	REJECTED,
	ACCEPTED,
	FAILED,
	FINISHED
};

/// @brief Represents a set of filter scores that have been accepted
class AcceptedScores : public utility::vector1< core::Real > {
public:
	AcceptedScores( core::Size const iteration ):
		utility::vector1< core::Real >(),
		iteration_( iteration ),
		rank_score_( 0.0 )
	{}
	AcceptedScores( core::Size const iteration, core::Real const rank_score ):
		utility::vector1< core::Real >(),
		iteration_( iteration ),
		rank_score_( rank_score )
	{}
	AcceptedScores( core::Size const iteration, core::Real const rank_score, utility::vector1< core::Real > const & scores ):
		utility::vector1< core::Real >( scores ),
		iteration_( iteration ),
		rank_score_( rank_score )
	{}
	friend std::ostream & operator<<( std::ostream & os, AcceptedScores const & scores );
	core::Size iteration() const { return iteration_; }
	core::Real rank_score() const { return rank_score_; }
	void set_rank_score( core::Real const rank_score ) { rank_score_ = rank_score; }
private:
	core::Size iteration_;
	core::Real rank_score_;
};

/// @brief GenericSimulatedAnnealer mover for performing simulated annealing trajectories
/// @details  This class extends GenericMonteCarloMover to provide temperature scaling
class GenericSimulatedAnnealer : public protocols::monte_carlo::GenericMonteCarloMover {
public:
	/// @brief default constructor
	GenericSimulatedAnnealer();

	/// @brief destructor
	~GenericSimulatedAnnealer() override;

	/// @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief apply GenericSimulatedAnnealer (Mover)
	void apply( Pose & pose ) override;

	/// @brief XML interface to this mover
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & ) override;

	/// @brief Clears accepted data and initialize with new pose
	virtual void reset( Pose & pose );

	/// @brief given a pose, score the result
	AcceptedScores
	score_pose( core::pose::Pose const & pose ) const;

public: // accessor
	/// @brief Saves the current state of the mover into checkpoint files
	void save_checkpoint_file() const;
	/// @brief Deletes the checkpoint files
	void remove_checkpoint_file() const;
	/// @brief tests to see if the checkpoint files exist and have been generated
	bool checkpoint_exists() const;
	/// @brief gets a list of scores for acceptance number i
	AcceptedScores const & accepted_scores( core::Size const i ) const { return accepted_scores_[i]; }
	/// @brief number of acceptances
	core::Size num_accepted_scores() const { return accepted_scores_.size(); }
	/// @brief if boltz_rank is used, this will calculate the ranking score from a list of filter scores
	core::Real calc_boltz_score( utility::vector1< core::Real > const & scores ) const;
	/// @brief create a silent struct tag for checkpointing with the given suffix
	std::string create_tag( std::string const & suffix ) const;

public: // mutators
	/// @brief loads checkpoint data and attempts to resume a run with pose as the pose
	void load_checkpoint_file( core::pose::Pose & pose );
	/// @brief set the checkpoint file name
	void checkpoint_file( std::string const & check_file ) { checkpoint_file_ = check_file; };
	/// @brief sets whether or not to delete checkpoint files when done -- useful for debugging
	void keep_checkpoint_file( bool const keep_checkpoint ) { keep_checkpoint_file_ = keep_checkpoint; };
	/// @brief scales temperatures by the factor provided
	void scale_temperatures( core::Real const temp_factor );
	/// @brief given a modified pose, determines whether we should accept or not, and updates internal class data accordingly
	///        uses randomly generated numbers to assess acceptance of scores with temperatures
	TrialResult boltzmann_result( core::pose::Pose & pose );
	/// @brief given a modified pose, determines whether we should accept or not, and updates internal class data accordingly
	TrialResult boltzmann_result( core::pose::Pose & pose,
		utility::vector1< core::Real > const & random_nums );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // private functions
	/// @brief calls a round of monte carlo -- basically copied from GenericMonteCarloMover
	TrialResult
	apply_mover( core::pose::Pose & pose );

	/// @brief computes the appropriately scaled temperatures and sets them
	void
	calculate_temps();

	utility::vector1< core::Real >
	calculate_standardized_scores( core::Size const filterid ) const;

	/// @brief calculates multiplier for temperatures based on which anneal step we are no
	core::Real calc_temp_factor() const;

	AcceptedScores
	read_checkpoint_line( std::istream & is ) const;

	void
	recompute_rank_scores();

private:
	/// @brief counter for how many accepts it takes to reach equilibrium at a given temperature
	core::Size history_;
	/// @brief mover to be called every eval_period steps
	protocols::moves::MoverOP periodic_mover_;
	/// @brief how many steps in between running the periodic mover?
	core::Size eval_period_;
	/// @brief what name should be used to store checkpoint information so that the run can be resumed?
	std::string checkpoint_file_;
	/// @brief if true, the checkpoint files will not be cleaned up after the apply() terminates. Useful for debugging/gathering statistics
	bool keep_checkpoint_file_;
	/// @brief a history of accepted scores, used for adapting temperature
	utility::vector1< AcceptedScores > accepted_scores_;
	/// @brief the initial temperatures, used for adapting temperature
	utility::vector1< core::Real > start_temperatures_;
	/// @brief counter that tells which annealing step we are on
	core::Size anneal_step_;
	/// @brief which step in this temperature are we on
	core::Size temp_step_;
	/// @brief what is the current trial number
	core::Size current_trial_;
};

/// @brief helper function that safely replaces a file with another
void replace_file( std::string const & origfile, std::string const & newfile );

} // namespace monte_carlo
} // namespace protocols

#endif

