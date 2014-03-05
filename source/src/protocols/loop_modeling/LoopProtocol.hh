// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopProtocol_HH
#define INCLUDED_protocols_loop_modeling_LoopProtocol_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>
#include <protocols/loop_modeling/loggers/Logger.fwd.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/loops/Loops.hh>

// Utility headers
#include <boost/utility.hpp>

namespace protocols {
namespace loop_modeling {

/// @brief Monte Carlo search for low energy loop conformations.
///
/// @details This class provides an easy way to run a Monte Carlo simulation 
/// searching for the lowest energy conformations for a set of loops.  This is 
/// most common way to use the classes in this namespace.  This simulation is 
/// organized into three nested loops.  On each iteration of the outermost 
/// loop, the lowest scoring pose is recovered and the repulsive terms in the 
/// score function may be ramped (although this ramping is disabled by 
/// default).  On each iteration of the intermediate loop, the temperature may 
/// be ramped (this ramping is enabled by default).  And on each iteration of 
/// the innermost loop, a new conformation is sampled and either accepted or 
/// rejected according to the Metropolis criterion.  The intermediate loop 
/// usually goes through more than 100 iterations, while the innermost and 
/// outermost loops only go through less than 5.
///
/// Like any mover, all the work is done by the apply() method.  The rest of 
/// the methods of this class are just getters and setters that can be used to 
/// control various aspects of the simulation.  The add_mover(), add_filter(), 
/// and add_acceptance_check() methods are worth drawing some attention to.  
/// These methods are used to build up the group of LoopMover objects that 
/// samples new loop conformations in the innermost loop.  The movers are 
/// guaranteed to be applied in the order they are added to the protocol.

class LoopProtocol : public LoopMover {

public:

	/// @brief Default constructor.
	LoopProtocol();

	/// @brief Default destructor.
	~LoopProtocol();

	/// @brief Return the class name of this mover.
	string get_name() const { return "LoopProtocol"; }

protected:

	/// @brief Use a Monte Carlo simulation to search for the best 
	/// conformations for the given loops.
	bool do_apply(Pose & pose);

private:

	/// @brief Setup the aspects of the simulation like the loop movers, the 
	/// score function, the fold tree, and the temperature.
	void start_protocol(Pose & pose);

	/// @brief Change the temperature of the simulation.  This is called once 
	/// every level-1 cycle.
	void ramp_temperature(Size iteration);

	/// @brief Change the weights of the repulsive score function terms.  This 
	/// is called once every level-2 cycle.
	void ramp_score_function(Size iteration);

	/// @brief Sample a new conformation by applying all the loop movers.  
	/// Afterwards apply the Metropolis accept-or-reject criterion.  This is 
	/// called once every level-3 cycle.
	void attempt_loop_move(Pose & pose, Size i, Size j, Size k);

	/// @brief Finalize any loggers and restore the original fold tree.
	void finish_protocol(Pose & pose);

public:

	/// @brief Add a LoopMover to the simulation.
	/// @details Loop movers will be applied in the order they're added.
	void add_mover(LoopMoverOP mover);

	/// @brief Add a Filter to the simulation.
	void add_filter(protocols::filters::FilterOP filter);

	/// @brief Add an acceptance check to the simulation.
	/// @details An acceptance check is always performed after all the loop 
	/// movers have been applied, unless one of the movers failed.  This method 
	/// allows additional acceptance checks to be included in the protocol by 
	/// added a loop mover that does nothing but make that check.
	void add_acceptance_check(string name="loop_move");

	/// @brief Add a logger to this simulation.
	/// @details The Logger subclasses are only meant to be used with 
	/// LoopProtocol.  These classes do not derive from LoopMover, and 
	/// implement a number of methods specifically designed for expressive 
	/// logging output.
	void add_logger(loggers::LoggerOP logger);

	/// @brief Set the number of iterations to use in this simulation.
	/// @details The simulation consists of three nested loops, so the given 
	/// list must contain three values.
	void set_iterations(IndexList values);

	/// @brief Set the number of iterations to use in this simulation.
	void set_iterations(Size i, Size j, Size k);

	/// @brief Set the initial and final temperatures for the simulation.
	/// @details The temperature will be linearly interpolated between these 
	/// values during the simulation.
	void set_temperature_schedule(Real initial, Real final);

	/// @brief Enable or disable temperature ramping in this simulation.
	void set_temperature_ramping(bool value);

	/// @brief Enable or disable score function ramping in this simulation.
	void set_score_function_ramping(bool value);

private:
	utilities::LoopMoverGroupOP movers_;
	protocols::moves::MonteCarloOP monte_carlo_;
	core::kinematics::FoldTree original_tree_;
	IndexList iterations_;
	loggers::LoggerList loggers_;

	bool ramp_score_function_;
	bool ramp_temperature_;

	Real initial_temperature_;
	Real final_temperature_;
	Real temperature_scale_factor_;
};

}
}

#endif
