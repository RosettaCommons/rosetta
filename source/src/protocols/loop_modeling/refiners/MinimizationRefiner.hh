// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_MinimizationRefiner_HH
#define INCLUDED_protocols_loop_modeling_refiners_MinimizationRefiner_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/simple_moves/MinMover.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

/// @brief Refine sampled loops using gradient minimization.
class MinimizationRefiner : public LoopMover {

public:

	/// @brief Constructor with options to configure the MinMover that will be
	/// used under the hood.
	MinimizationRefiner(
		bool cartesian=false,
		core::optimization::MinimizerOptionsOP options=
		core::optimization::MinimizerOptionsOP());

	/// @copydoc LoopMover::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "MinimizationRefiner"; }

	/// @brief Specify whether or not Cartesian minimization should be used.
	void use_cartesian(bool setting);

	/// @brief Return true if cartesian minimization will be used.  The
	/// alternative is atom tree minimization.
	bool use_cartesian() const;

	/// @brief Get the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP score_function);

	/// @brief Set the minimizer options.
	void set_min_options(core::optimization::MinimizerOptionsOP options);

	/// @brief Non-const access to the minimizer options.  May be NULL.
	core::optimization::MinimizerOptionsOP get_min_options();

	/// @brief Const access to the minimizer options.  May be NULL.
	core::optimization::MinimizerOptionsCOP get_min_options() const;

protected:

	/// @brief Perform Cartesian minimization within 10A of the loops being
	/// sampled.
	bool do_apply(Pose & pose);

private:

	bool use_cartesian_;
	core::optimization::MinimizerOptionsOP min_options_;

};

}
}
}

#endif

