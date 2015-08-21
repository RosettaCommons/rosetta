// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh
/// @brief Abstract class to define interface for all types of "inner cycle" operations used for loop refinement.
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycle_HH
#define INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycle_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
// #include <protocols/loops/loop_mover/LoopMover.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class LoopRefineInnerCycle: public moves::Mover {
public: // boiler plate / virtuals
	// default constructor
	LoopRefineInnerCycle();

	// copy constructor
	LoopRefineInnerCycle( LoopRefineInnerCycle const & rhs );

	// assignment operator
	LoopRefineInnerCycle & operator=( LoopRefineInnerCycle const & rhs );

	// destructor
	virtual ~LoopRefineInnerCycle();

	// constructor with arguments
	LoopRefineInnerCycle(
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	);

	virtual void apply( Pose & ) = 0;
	virtual std::string get_name() const;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	/// @brief Associates relevant options with the LoopRefineInnerCycle class
	static void register_options();

	// NOTE: The clone() and fresh_instance() virtual methods are omitted because this class is abstract

public: // printing methods
	virtual void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, LoopRefineInnerCycle const & loop_refine_inner_cycle );

public: // class-specific public methods
	bool debug() const;
	void set_debug( bool debug );

	moves::MonteCarloOP mc() const;
	virtual void set_mc( moves::MonteCarloOP mc);

	core::scoring::ScoreFunctionOP scorefxn() const;
	virtual void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn);

	core::pack::task::TaskFactoryOP task_factory() const;
	virtual void set_task_factory( core::pack::task::TaskFactoryOP tf );

	LoopMover_Refine_CCDAP loop_mover() const;
	virtual void set_loop_mover( LoopMover_Refine_CCDAP new_owner_in_town );

protected:
	core::kinematics::MoveMapOP movemap() const;
	void set_movemap( core::kinematics::MoveMapOP movemap );

	Loops get_one_random_loop() const;

	void setup_objects( Pose const & pose );

private: // methods
	void init();
	void init(
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	);

	void init_for_equal_operator_and_copy_constructor( LoopRefineInnerCycle & lhs, LoopRefineInnerCycle const & rhs);
	void init_options();

private: // data
	bool debug_;

	LoopMover_Refine_CCDAP loop_mover_that_owns_me_; // access pointer prevents strong reference cycle
	moves::MonteCarloOP mc_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pack::task::TaskFactoryOP tf_;

	mutable core::kinematics::MoveMapOP movemap_;

}; // class LoopRefineInnerCycle

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycle_HH
