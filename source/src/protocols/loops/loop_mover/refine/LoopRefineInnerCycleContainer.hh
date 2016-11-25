// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainer.hh
/// @brief This class is a LoopRefineInnerCycle that contains one or more other LoopRefineInnerCycles to allow a developer to
/// quickly string together existing LoopRefineInnerCycles in new ways to create new loop refinement protocols.
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycleContainer_HH
#define INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycleContainer_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleContainer.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.hh>

// Package headers

// Project headers
#include <protocols/moves/MoverContainer.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class LoopRefineInnerCycleContainer: public LoopRefineInnerCycle {
public: // boiler plate / virtuals
	// default constructor
	LoopRefineInnerCycleContainer();

	// copy constructor
	LoopRefineInnerCycleContainer( LoopRefineInnerCycleContainer const & rhs );

	// assignment operator
	LoopRefineInnerCycleContainer & operator=( LoopRefineInnerCycleContainer const & rhs );

	// destructor
	virtual ~LoopRefineInnerCycleContainer();

	void apply( Pose & ) override;
	// XRW TEMP  virtual std::string get_name() const;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	bool reinitialize_for_new_input() const override;

	/// @brief Associates relevant options with the LoopRefineInnerCycleContainer class
	static void register_options();

public: // printing methods
	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, LoopRefineInnerCycleContainer const & loop_refine_inner_cycle_container );

public: // class-specific public methods
	void add_inner_cycle_step( LoopRefineInnerCycleOP inner_cycle_step );

	// overriden methods: Not 'virtual' in the sense that subclasses ought to override them. This is a special case.
	void set_mc( moves::MonteCarloOP mc) override;
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn) override;
	void set_task_factory( core::pack::task::TaskFactoryOP tf ) override;
	void set_loop_mover( LoopMover_Refine_CCDAP new_owner_in_town ) override;

	// This one comes from moves::Mover
	void set_native_pose( PoseCOP pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods
	void setup_objects( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( LoopRefineInnerCycleContainer & lhs, LoopRefineInnerCycleContainer const & rhs);
	void init_options();

private: // data
	typedef utility::vector1< LoopRefineInnerCycleOP > InnerCycleList;
	InnerCycleList inner_cycle_list_;
	moves::SequenceMoverOP inner_cycle_steps_;

}; // class LoopRefineInnerCycleContainer

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_mover_refine_LoopRefineInnerCycleContainer_HH
