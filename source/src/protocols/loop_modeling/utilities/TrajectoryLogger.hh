// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_TrajectoryLogger_HH
#define INCLUDED_protocols_loop_modeling_utilities_TrajectoryLogger_HH

// Unit headers
#include <protocols/loop_modeling/utilities/TrajectoryLogger.fwd.hh>
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopBuilder.fwd.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <basic/Tracer.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Log a detailed account of everything that happens in a loop modeling
/// simulation.
/// @details Basic usage is to initialize the logger either with one of the
/// constructors or the init() method, then to call the record methods (e.g.
/// record_move()) whenever something of note happens.  The record method were
/// written with LoopProtocol and LoopBuilder in mind, so you may have to add
/// some of your own record methods if you want to use this logger with another
/// LoopMover subclass.
class TrajectoryLogger {

public:

	/// @brief Construct an empty logger.
	/// @details Use the init() method to give the setup the logger before any of
	/// the record methods are called.
	TrajectoryLogger(
		string const & prefix="");

	/// @brief Construct a logger from a loop mover.
	/// @details The given loop mover must have both a loops and a score function
	/// at the time to constructor is called.  If either condition is not met, an
	/// exception will be raised.
	TrajectoryLogger(
		LoopMover const * mover,
		string const & prefix="");

	/// @brief Construct a logger from a score function and a loops object.
	TrajectoryLogger(
		protocols::loops::LoopsCOP loops,
		core::scoring::ScoreFunctionCOP scorefxn,
		string const & prefix="");

	/// @brief Initialize the logger from a loop mover after it has been
	/// constructed.
	void init(
		LoopMover const * mover);

	/// @brief Initialize the logger from a score function and a loops object
	/// after it has been constructed.
	void init(
		protocols::loops::LoopsCOP loops,
		core::scoring::ScoreFunctionCOP scorefxn);

	/// @brief Record a single Monte Carlo move.
	/// @details This overload is meant for use in LoopProtocol, which uses three
	/// nested for-loops to ramp some score function terms and the temperature at
	/// different rates.
	/// @param tr The output stream to record the move to.
	/// @param pose The pose on which the last move was made (used for
	/// calculating a score and an RMSD).
	/// @param i The current "score function" cycle.
	/// @param j The current "temperature" cycle.
	/// @param k The current "mover" cycle.
	/// @param proposed Whether or not a move was proposed (i.e. did the loop
	/// mover succeed in producing a new conformation).
	/// @param accepted Whether or not the last move was accepted (i.e. did it
	/// pass the Metropolis criterion).
	void record_move(
		basic::Tracer & tr,
		Pose & pose,
		Size const i,
		Size const j,
		Size const k,
		bool const proposed,
		bool const accepted) const;

	/// @brief Record a single Monte Carlo move.
	/// @details This overload is meant for use in LoopBuilder, which generates
	/// loop conformation is a single non-nested for-loop.
	/// @param tr The output stream to record the move to.
	/// @param pose The pose on which the last move was made (used for
	/// calculating a score and an RMSD).
	/// @param i The current iteration.
	/// @param proposed Whether or not a move was proposed (i.e. did the loop
	/// mover succeed in producing a new conformation).
	/// @param accepted Whether or not the last move was accepted (i.e. did it
	/// pass the Metropolis criterion).
	void record_move(
		basic::Tracer & tr,
		Pose & pose,
		Size const iteration,
		bool const proposed,
		bool const accepted) const;

	/// @brief Record when a low-scoring pose is recovered.
	void record_new_pose(
		basic::Tracer & tr,
		Pose & pose) const;

	/// @brief Record when a score function term is ramped (and scores change as
	/// a result).
	void record_new_score_function(
		basic::Tracer & tr,
		Pose & pose) const;

	/// @brief Record when the temperature is ramped.
	void record_new_temperature(
		basic::Tracer & tr,
		Real temperature) const;

	/// @brief Record the end of a protocol.
	/// @details Print out some more expensive quality metrics and attach some
	/// scores to the pose itself.
	void record_endpoint(
		basic::Tracer & tr,
		Pose & pose) const;

private:

	/// @brief Return the score of the given pose.
	Real calc_score(Pose & pose) const;

	/// @brief Return the heavy-atom backbone RMSD between the given pose and the
	/// <tt>-in:%file:native</tt> pose within the given loops.
	Real calc_rmsd_to_native(Pose & pose) const;

public:

	/// @brief Get the prefix that will be used for the scores associated with
	/// the pose at the end of the trajectory.
	string get_prefix() const;

	/// @brief Set the prefix that will be used for the scores associated with
	/// the pose at the end of the trajectory.
	void set_prefix(string);

	/// @brief Get the reference structure for RMSD calculations.
	/// @see native_pose_
	Pose const & get_native_pose() const;

	/// @brief Set the reference structure for RMSD calculations.
	/// @see native_pose_
	void set_native_pose(Pose const &);

	/// @brief Disable RMSD calculations.
	void unset_native_pose();

	/// @brief Get the region of the protein where an backbone heavy-atom RMSD
	/// will be calculated for each step in the trajectory.
	/// @see loops_
	protocols::loops::LoopsCOP get_loops() const;

	/// @brief Set the region of the protein where an backbone heavy-atom RMSD
	/// will be calculated for each step in the trajectory.
	/// @see loops_
	void set_loops(protocols::loops::LoopsCOP);

	/// @brief Get the score function used to score each step in the trajectory.
	core::scoring::ScoreFunctionCOP get_score_function() const;

	/// @brief Set the score function used to score each step in the trajectory.
	void set_score_function(core::scoring::ScoreFunctionCOP);

	/// @brief Return how long the trajectory has been running in seconds.
	long get_timer() const;

	/// @brief Indicate that the trajectory is starting now.
	/// @details The timer is automatically started when the logger is
	/// constructed, you only need to call this method if the logger is
	/// constructed significantly before the trajectory in question begins.
	void reset_timer();

private:

	/// @brief The prefix that will be used for the scores associated with the
	/// pose at the end of the trajectory.
	string prefix_;

	/// @brief The reference structure for RMSD calculations.
	/// @details This structure is typically initialized from
	/// <tt>-in:%file:native</tt> in the constructor.  Reading a pose in from a
	/// file is fairly expensive, and we want to calculate RMSDs after every
	/// iteration, so it's important that this initialization only happen once.
	Pose native_pose_;

	/// @brief Indicate that RMSDs can be calculated because there is a native
	/// pose.
	bool have_native_pose_ = false;

	/// @brief The region of the protein where an backbone heavy-atom RMSD will
	/// be calculated for each step in the trajectory.
	/// @details The RMSD is calculated relative to the <tt>-in:%file:native</tt>
	/// structure.  If that option was not specified, no RMSD will be calculated.
	protocols::loops::LoopsCOP loops_ = nullptr;

	/// @brief The score function used to score each step in the trajectory.
	core::scoring::ScoreFunctionCOP scorefxn_ = nullptr;

	/// @brief The number of seconds since the epoch that had elapsed at the
	/// beginning of the trajectory.
	long start_time_ = 0;

};

} // namespace utilities
} // namespace loop_modeling
} // namespace protocols

#endif
