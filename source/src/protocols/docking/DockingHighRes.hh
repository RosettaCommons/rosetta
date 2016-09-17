// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file docking_initialization_protocols
/// @brief initialization protocols for docking
/// @details
///  This contains the functions that create initial positions for docking
///  You can either randomize partner 1 or partner 2, spin partner 2, or
///  perform a simple perturbation.
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov


#ifndef INCLUDED_protocols_docking_DockingHighRes_hh
#define INCLUDED_protocols_docking_DockingHighRes_hh

#include <protocols/docking/DockingHighRes.fwd.hh>

// Package headers
#include <protocols/docking/types.hh>
#include <protocols/docking/DockTaskFactory.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>


#include <utility/vector1.hh>


// option key includes

namespace protocols {
namespace docking {

class DockingHighRes : public moves::Mover
{
	typedef core::Real Real;
public:

	/// @brief Default constructor
	DockingHighRes();

	/// @brief Constructor with one argument - the jump number to dock over.
	DockingHighRes(
		core::Size const rb_jump
	);

	/// @brief Constructor with two arguments.  The first is thejump number, the second is a scorefunction that will be
	///  used for docking and packing.
	DockingHighRes(
		core::Size const rb_jump,
		core::scoring::ScoreFunctionOP scorefxn
	);

	/// @brief Constructor with three arguments.  The first is the jump number, the second is a scorefunction that will
	///  be used for docking and the third is a scorefunction that will be used for packing.
	DockingHighRes(
		core::Size const rb_jump,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	/// @brief Constructor with three arguments.  The first is the DockJumps, the second is a scorefunction that will
	///  be used for docking and the third is a scorefunction that will be used for packing.
	DockingHighRes(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	//Copy constructor
	DockingHighRes( DockingHighRes const & old_instance );

	~DockingHighRes() override;

	//clone
	//protocols::moves::MoverOP clone() const = 0;

	void set_task_factory( core::pack::task::TaskFactoryCOP tf );
	//allows one to ignore the DockTaskFactoryOP and allow use of whatever is given in set_task_factory()
	//only works for DockMCMProtocol at the moment
	void set_ignore_default_task( bool ignore_default_task );
	bool ignore_default_task();
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void set_scorefxn_pack( core::scoring::ScoreFunctionOP scorefxn_pack );

	DockJumps & movable_jumps() { return movable_jumps_; }
	DockJumps const & movable_jumps() const { return movable_jumps_; }
	void set_movable_jumps( DockJumps movable_jumps ) { movable_jumps_ = movable_jumps; }
	core::pack::task::TaskFactoryCOP task_factory();

	// protocol functions
	void apply( core::pose::Pose & pose ) override = 0;
	std::string get_name() const override = 0;

	void set_sc_min( bool sc_min ){ sc_min_ = sc_min; }
	void set_rt_min( bool rt_min ){ rt_min_ = rt_min; }
	void set_partners( std::string partners ) { partners_ = partners; }
  std::string get_partners( ) const { return partners_; }
	void set_interface_definition_task_operation( protocols::toolbox::task_operations::InterfaceTaskOperationOP interface_definition );
	void set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations );
	void add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation );
	utility::vector1< core::pack::task::operation::TaskOperationOP > get_additional_task_operarations();


	bool sc_min() { return sc_min_; }
	bool rt_min() { return rt_min_; }
	std::string partners() { return partners_; }

	core::scoring::ScoreFunctionOP scorefxn() const;
	core::scoring::ScoreFunctionOP scorefxn_pack() const;

protected:
	protocols::docking::DockTaskFactoryOP tf2();  //JQX: change COP to OP

private:
	bool sc_min_;
	bool rt_min_;
	std::string partners_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_pack_;

	// docking
	DockJumps movable_jumps_;

	// the task factory that will be used for all parts of docking
	core::pack::task::TaskFactoryOP tf_;
	protocols::docking::DockTaskFactoryOP tf2_; // tf2 = task factory factory - we use this to create our task factories
	bool ignore_default_task_; //will ignore the DockTaskFactoryOP and allow use of whatever is given in set_task_factory()
	void init( DockJumps const movable_jumps );
};

} // docking
} // protocols

#endif
