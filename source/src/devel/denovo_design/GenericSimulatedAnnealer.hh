// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/devel/simple_moves/GenericSimulatedAnnealer.hh
/// @brief perform a given mover and sample structures by MonteCarlo with simulated annealing of acceptance temperatures
/// @detailed The "score" evaluation of pose during MC after applying mover is done by
/// ither FilterOP that can do report_sm() or ScoreFunctionOP you gave.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_devel_denovo_design_GenericSimulatedAnnealer_hh
#define INCLUDED_devel_denovo_design_GenericSimulatedAnnealer_hh

// Unit header
#include <devel/denovo_design/GenericSimulatedAnnealer.fwd.hh>
#include <protocols/moves/DataMapObj.hh>

// Package headers
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/simple_moves/GenericMonteCarloMover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// External headers

namespace devel {
namespace denovo_design {

enum TrialResult {
		REJECTED,
		ACCEPTED,
		FAILED,
		FINISHED
};

class GenericSimulatedAnnealer : public protocols::simple_moves::GenericMonteCarloMover {
public:
	/// @brief default constructor
	GenericSimulatedAnnealer();

	/// @brief destructor
	~GenericSimulatedAnnealer();

	/// @brief create copy constructor
	virtual protocols::moves::MoverOP clone() const;

	/// @brief create this type of objectt
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief apply GenericSimulatedAnnealer (Mover)
	virtual void apply( Pose & pose );

	virtual std::string get_name() const;

	virtual void parse_my_tag(
		TagPtr const tag,
		DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & );

	virtual void reset( Pose & pose );

public: // accessor

public: // mutators

private: // private functions
	/// @brief calls a round of monte carlo -- basically copied from GenericMonteCarloMover
	TrialResult
	apply_mover( core::pose::Pose & pose );

	/// @brief computes the appropriately scaled temperatures and sets them
	void
	calculate_temps();

	utility::vector1< core::Real >
	calculate_standardized_scores( core::Size const filterid ) const;

	/// @brief scales temperatures by the factor provided
	void
	scale_temperatures( core::Real const temp_factor );

private:
	/// @brief counter for how many accepts it takes to reach equilibrium at a given temperature
	core::Size history_;
	/// @brief a history of accepted scores, used for adapting temperature
	utility::vector1< utility::vector1< core::Real > > accepted_scores_;
	/// @brief the initial temperatures, used for adapting temperature
	utility::vector1< core::Real > start_temperatures_;
	/// @brief counter that tells which annealing step we are on
	core::Size anneal_step_;
	core::Size temp_step_;
};

} // namespace denovo_design
} // namespace devel

#endif

