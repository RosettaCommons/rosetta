// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopProtocol_HH
#define INCLUDED_protocols_loop_modeling_LoopProtocol_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.fwd.hh>
#include <protocols/loop_modeling/LoopModelerTests.fwd.hh> //For friendship

// Core headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocols headers
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/loops/Loops.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <boost/utility.hpp>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace loop_modeling {

/// @brief Monte Carlo search for low energy loop conformations.
///
/// @details This class provides an easy way to run a Monte Carlo simulation
/// searching for the lowest energy conformations for a set of loops.  This
/// simulation is organized into three nested loops.  On each iteration of the
/// outermost loop, the lowest scoring pose is recovered and the repulsive
/// terms in the score function may be ramped (although this ramping is
/// disabled by default).  On each iteration of the intermediate loop, the
/// temperature may be ramped (this ramping is enabled by default).  And on
/// each iteration of the innermost loop, a new conformation is sampled and
/// either accepted or rejected according to the Metropolis criterion.  The
/// intermediate loop usually goes through more than 100 iterations, while the
/// innermost and outermost loops only go through less than 5.
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
	~LoopProtocol() override;

	/// @copydoc LoopMover::get_name
	// XRW TEMP  string get_name() const override { return "LoopProtocol"; }

protected:

	/// @brief Use a Monte Carlo simulation to search for the best
	/// conformations for the given loops.
	bool do_apply(Pose & pose) override;

public:

	/// @copydoc LoopMover::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose) override;

private:

	/// @brief Setup the aspects of the simulation like the loop movers, the
	/// score function, the fold tree, and the temperature.
	void start_protocol(Pose & pose);

	/// @brief Change the weights of the repulsive score function terms.  This is
	/// called once every level-2 cycle.
	void ramp_score_function(Size iteration);

	/// @brief Change the temperature of the simulation.  This is called once
	/// every level-1 cycle.
	void ramp_temperature(Size iteration);

	/// @brief Sample a new conformation by applying all the loop movers.
	/// Afterwards apply the Metropolis accept-or-reject criterion.  This is
	/// called once every level-3 cycle.
	void attempt_loop_move(Pose & pose, Size i, Size j, Size k);

	/// @brief No-op right now, but may be useful in the future.
	void finish_protocol(Pose & pose);

public:

	/// @brief Add a LoopMover to the simulation.
	/// @details Loop movers will be applied in the order they're added.
	/// @see add_refiner()
	/// @see mark_as_default()
	void add_mover(LoopMoverOP mover);

	/// @brief Add a LoopMover to the simulation as a refiner.
	/// @details This method is very similar to add_mover().  Both methods add a
	/// new LoopMover to the simulation.  This difference is subtle, but relevant
	/// when a LoopProtocol is being filled with a default set of movers that
	/// might be modified later.
	///
	/// This situation arises when a LoopModeler object is being constructed, or
	/// when a <LoopModeler> tag is being parsed by rosetta scripts.  In both of
	/// these cases, centroid and fullatom LoopProtocols are created right away
	/// and filled with default movers.  This avoids the necessity of specifying
	/// all the usual movers every time loop modeling is invoked.  But the
	/// default movers may be subsequently overridden.  In this case, there is a
	/// difference between LoopMovers added via add_mover() (i.e. movers) and
	/// add_refiner() (i.e. refiners).
	///
	/// The difference is that movers can be marked as "default" and then be
	/// implicitly overridden.  Consider the following sequence of calls:
	///
	///     protocol->add_mover(new CcdMover);
	///     protocol->mark_as_default();
	///     protocol->add_mover(new KicMover);
	///
	/// This will result in a LoopProtocol containing only a KicMover.  When
	/// mark_as_default() is called, all movers within the protocol are set to be
	/// replaced the next time add_mover() is called.  You can think of those
	/// movers as being part of a default configuration that should be replaced
	/// when a new configuration is given.
	///
	/// Refiners are not affected by mark_as_default().  Calling add_refiner()
	/// will always add a refiner and will never remove old ones, although old
	/// refiners can be explicitly removed with clear_refiners().  Refiners are
	/// also applied after movers, regardless of the order in which add_mover()
	/// and add_refiner() were called.  The purpose of all this is to make it
	/// easy to replace the parts of a LoopProtocol that are more variable
	/// without having to also replace those that are more static.
	///
	/// The variable part of a LoopProtocol is typically the sampling algorithm:
	/// KIC or CCD, often with myriad options.  The static part is typically the
	/// refinement algorithm: minimization, rotamer trials, and repacking, all
	/// with carefully tuned parameters.  Using add_mover(), mark_as_default(),
	/// and add_refiner() makes it easy to change the sampling algorithm without
	/// having to worry about the refinement algorithm.
	/// @see add_mover()
	/// @see mark_as_default()
	void add_refiner(LoopMoverOP refiner);

	/// @brief Add a Filter to the simulation.
	void add_filter(protocols::filters::FilterOP filter);

	/// @brief Add an acceptance check to the simulation.
	/// @details An acceptance check is always performed after all the loop
	/// movers have been applied, unless one of the movers failed.  This method
	/// allows additional acceptance checks to be included in the protocol by
	/// adding a loop mover that does nothing but make that check.
	void add_acceptance_check(string name="loop_move");

	/// @brief Remove all movers from the protocol.
	/// @see add_mover()
	/// @see add_refiner()
	/// @see clear_refiners()
	/// @see clear_movers_and_refiners()
	void clear_movers();

	/// @brief Remove all refiners from the protocol.
	/// @see add_mover()
	/// @see add_refiner()
	/// @see clear_movers()
	/// @see clear_movers_and_refiners()
	void clear_refiners();

	/// @brief Remove all movers and refiners from the protocol.
	/// @see add_mover()
	/// @see add_refiner()
	/// @see clear_movers()
	/// @see clear_refiners()
	void clear_movers_and_refiners();

	/// @brief Indicate that every mover that is currently part of the protocol
	/// should be removed the next time a mover is added.
	/// @see add_mover()
	/// @see add_refiner()
	void mark_as_default();

	/// @brief Get the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP score_function);

	/// @brief Set the number of iterations the loop that ramps the score
	/// function should make.
	/// @details The score function loop is the outermost loop.  It contains the
	/// temperature loop and the mover loop.
	void set_sfxn_cycles(Size x);

	/// @brief Get the number of iterations the loop that ramps the score
	/// function will make.
	Size get_sfxn_cycles();

	/// @brief Set the number of iterations the loop that ramps the temperature
	/// should make.
	/// @details The temperature loop is contained by the score function loop.
	/// It contains the mover loop.  If the second boolean argument is true, the
	/// number of iterations will be the first value times the combined length of
	/// all the loops being sampled.
	void set_temp_cycles(Size x, bool times_loop_length=false);

	/// @brief Get the number of iterations the loop that ramps the temperature
	/// should make.
	/// @details This will throw an exception if the 'times_loop_length' argument
	/// was passed to set_temp_cycles() and no loop has been specified yet.  The
	/// value returned will always be the actual number of iterations that will
	/// be carried out.
	Size get_temp_cycles();

	/// @brief Specify how many times the loops movers should be invoked after
	/// the score function and temperature have been updated.
	/// @details The mover loop is the innermost loop.  It is contained by the
	/// score function loop and the temperature loop.
	void set_mover_cycles(Size x);

	/// @brief Return how many times the loops movers will be invoked after the
	/// score function and temperature have been updated.
	Size get_mover_cycles();

	/// @brief Indicate that this protocol is being used for a test run.
	/// @details This method simply limits the number of iterations that will be
	/// performed.  It can be "undone" by explicitly setting the score function,
	/// temperature, and mover cycles back to normal values.
	void mark_as_test_run();

	/// @brief Enable or disable ramping the repulsive term of the score
	/// function.
	void set_repulsive_term_ramping(bool value);

	/// @brief Enable or disable ramping the repulsive term of the score
	/// function.
	void set_rama_term_ramping(bool value);

	/// @brief Enable or disable temperature ramping in this simulation.
	void set_temperature_ramping(bool value);

	/// @brief Set the initial and final temperatures for the simulation.
	/// @details The temperature will be linearly interpolated between these
	/// values during the simulation.
	void set_temperature_schedule(Real start, Real stop);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	friend class ::LoopModelerTests;

private:
	utilities::LoopMoverGroupOP protocol_;
	utilities::LoopMoverGroupOP movers_, refiners_;
	protocols::moves::MonteCarloOP monte_carlo_;

	Size sfxn_cycles_;
	Size temp_cycles_;
	Size mover_cycles_;

	bool ramp_sfxn_rep_;
	bool ramp_sfxn_rama_;
	bool ramp_temp_;
	bool scale_temp_cycles_;

	Real initial_temp_;
	Real final_temp_;
	Real original_rama_weight_;
	Real original_rama2b_weight_;
	Real original_repulsive_weight_;

	bool test_run_;
};

}
}

#endif
