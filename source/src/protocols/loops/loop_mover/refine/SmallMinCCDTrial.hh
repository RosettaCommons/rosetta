// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/refine/SmallMinCCDTrial.hh
/// @brief Perform a small move followed CCD closure, packing and minimization
/// @detailed
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_loops_loop_mover_refine_SmallMinCCDTrial_HH
#define INCLUDED_protocols_loops_loop_mover_refine_SmallMinCCDTrial_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/SmallMinCCDTrial.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>

// Package headers

// Project headers
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility headers

// C++ headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class SmallMinCCDTrial : public LoopRefineInnerCycle {
public: // boiler plate / virtuals
	// default constructor
	SmallMinCCDTrial();

	// copy constructor
	SmallMinCCDTrial( SmallMinCCDTrial const & rhs );

	// assignment operator
	SmallMinCCDTrial & operator=( SmallMinCCDTrial const & rhs );

	// destructor
	virtual ~SmallMinCCDTrial();

	// constructor with arguments
	SmallMinCCDTrial(
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	);

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	/// @brief Associates relevant options with the LoopRefineInnerCycle class
	static void register_options();
	
public: // printing methods
	virtual void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, SmallMinCCDTrial const & small_min_ccd_trial );

public: // class-specific public methods
	core::Size number_of_moves() const;
	void set_number_of_moves( core::Size nmoves );

	core::optimization::MinimizerOptionsOP minimizer_options() const;
	void set_minimizer_options( core::optimization::MinimizerOptionsOP minimizer_options );
	
private: // methods
	void setup_objects( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( SmallMinCCDTrial & lhs, SmallMinCCDTrial const & rhs);
	void init_options();

	core::optimization::AtomTreeMinimizerOP minimizer( core::pose::Pose const & pose ) const;

private: // data
	core::Size nmoves_;
	core::optimization::MinimizerOptionsOP minimizer_options_;
	mutable core::optimization::AtomTreeMinimizerOP minimizer_;

private: // Excessive debugging output
	void debug_zero( Pose & pose );
	void debug_one( Pose & pose );
	void debug_two( Pose & pose );
	void debug_three( Pose & pose );
	void debug_four( Pose & pose );
	void debug_five( Pose & pose );

}; // class SmallMinCCDTrial

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_SmallMinCCDTrial_HH
