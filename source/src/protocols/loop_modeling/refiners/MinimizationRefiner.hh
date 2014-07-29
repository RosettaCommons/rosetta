// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_MinimizationRefiner_HH
#define INCLUDED_protocols_loop_modeling_refiners_MinimizationRefiner_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.fwd.hh>

// Core headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.fwd.hh>

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
			core::optimization::MinimizerOptionsOP options=NULL);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "MinimizationRefiner"; }

	/// @brief Specify whether or not Cartesian minimization should be used.
	void use_cartesian(bool setting);

	/// @brief Return true if cartesian minimization will be used.  The 
	/// alternative is AtomTree minimization.
	bool use_cartesian() const;

	/// @brief Set the minimizer options.
	void min_options(core::optimization::MinimizerOptionsOP options);
	
	/// @brief Non-const access to the minimizer options.  May be NULL.
	core::optimization::MinimizerOptionsOP min_options();
	
	/// @brief Const access to the minimizer options.  May be NULL.
	core::optimization::MinimizerOptionsCOP min_options() const;
	
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

