// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/refine/ShearMinCCDTrial.hh
/// @brief Concrete class derived from LoopRefineInnerCycle to implement the CCD min trial flavor of inner cycle refinement.
/// @details
///
/// @author Michael Pacella (mpacella88@gmail.com)


#ifndef INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_HH
#define INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/ShearMinCCDTrial.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>

// Package headers
#include <protocols/loops/loop_mover/LoopMover.fwd.hh>

// Project headers
//#include <core/pack/task/TaskFactory.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>


// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class ShearMinCCDTrial: public LoopRefineInnerCycle {
public: // boiler plate / virtuals
	// default constructor
	ShearMinCCDTrial();

	// constructor with arguments
	ShearMinCCDTrial(
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	);

	// copy constructor
	ShearMinCCDTrial( ShearMinCCDTrial const & rhs );

	// assignment operator
	ShearMinCCDTrial & operator=( ShearMinCCDTrial const & rhs );

	// destructor
	~ShearMinCCDTrial();

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	// clone and fresh instance
	virtual moves::MoverOP clone() const;

	virtual moves::MoverOP fresh_instance() const;

public: // printing methods
	virtual void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, ShearMinCCDTrial const & loop_refine_shear_CCD_min_trial_inner_cycle );

public: // class-specific public methods

private: // methods
	void setup_objects( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( ShearMinCCDTrial & lhs, ShearMinCCDTrial const & rhs);
	void init_options();
	core::optimization::AtomTreeMinimizerOP minimizer( core::pose::Pose const & pose ) const;

private: //data
	Size nmoves_;
	core::optimization::MinimizerOptionsOP min_options_;
	mutable core::optimization::AtomTreeMinimizerOP minimizer_;

}; // class ShearMinCCDTrial

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_ShearMinCCDTrial_HH
