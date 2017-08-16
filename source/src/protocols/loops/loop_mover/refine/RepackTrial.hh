// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/RepackTrial.hh
/// @brief Abstract class to define interface for all types of "inner cycle" operations used for loop refinement.
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_loops_loop_mover_refine_RepackTrial_HH
#define INCLUDED_protocols_loops_loop_mover_refine_RepackTrial_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/RepackTrial.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>

// Package headers

// Project headers
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility headers

// C++ headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class RepackTrial : public LoopRefineInnerCycle {
public: // boiler plate / virtuals
	// default constructor
	RepackTrial();

	// copy constructor
	RepackTrial( RepackTrial const & rhs );

	// assignment operator
	RepackTrial & operator=( RepackTrial const & rhs );

	// destructor
	virtual ~RepackTrial();

	// constructor with arguments
	RepackTrial(
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	);

	void apply( Pose & ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	bool reinitialize_for_new_input() const override;

	/// @brief Associates relevant options with the LoopRefineInnerCycle class
	static void register_options();

public: // printing methods
	void show( std::ostream & out=std::cout ) const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	friend std::ostream & operator<<(std::ostream& out, RepackTrial const & repack_trial );

public: // class-specific public methods

private: // methods
	void setup_objects( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( RepackTrial & lhs, RepackTrial const & rhs);
	void init_options();

private: // data


}; // class RepackTrial

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_RepackTrial_HH
