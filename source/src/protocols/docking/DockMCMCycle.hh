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


#ifndef INCLUDED_protocols_docking_DockMCMCycle_hh
#define INCLUDED_protocols_docking_DockMCMCycle_hh

#include <protocols/docking/types.hh>
#include <protocols/docking/DockMCMCycle.fwd.hh>

// Package headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/TrialMover.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/MoverContainer.fwd.hh>
#include <utility/vector1.hh>
#include <iostream>


// option key includes

namespace protocols {
namespace docking {

class DockMCMCycle : public moves::Mover
{
	typedef core::Real Real;
public:

	/// @brief Default constructor
	DockMCMCycle();

	/// @brief Constructor with two arguments: The first argument is the jump number to dock over.
	///  The second is a scorefunction that will be used for docking and packing.
	DockMCMCycle(
		core::Size const rb_jump,
		core::scoring::ScoreFunctionOP scorefxn
	);

	/// @brief Constructor with arguments: The first argument is the jump number to dock over.
	///  The second is a scorefunction that will be used for docking, the third is a scorefunction that will be used for packing.
	DockMCMCycle(
		core::Size const rb_jump,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	/// @brief Constructor with arguments: The first argument is a DockJumps vector.
	///  The second is a scorefunction that will be used for docking, the third is a scorefunction that will be used for packing.
	DockMCMCycle(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	~DockMCMCycle() override;

	//clone
	protocols::moves::MoverOP clone() const override;

	/// @brief Sets the default values for all members
	void set_default();

	void set_move_map( core::kinematics::MoveMapOP movemap );
	void set_min_type( std::string min_type ) { min_type_ = min_type; }
	void set_task_factory( core::pack::task::TaskFactoryCOP tf );

	moves::MonteCarloOP get_mc() { return mc_; }

	// protocol functions
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, const DockMCMCycle & dp );


	//JQX: allow the DockMCMProtocol.cc file to change the index of the CycleMover
	void reset_cycle_index();
	void init_mc(core::pose::Pose & pose);

	void set_scmin(bool setting){scmin_=setting;}
	void set_rtmin(bool setting){rtmin_=setting;}
	void set_rot_magnitude(core::Real value){rot_magnitude_=value;}

	DockJumps const & movable_jumps() const;

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_pack_;
	core::kinematics::MoveMapOP movemap_;
	// moves::TrialMoverOP dock_mcm_mover_;
	moves::CycleMoverOP dock_mcm_cycle_; //JQX define it
	moves::MonteCarloOP mc_;


	// docking
	DockJumps movable_jumps_;
	core::Real trans_magnitude_, rot_magnitude_;
	bool rtmin_, scmin_; // belongs to the packer task

	void setup_protocol( core::pose::Pose & pose );
	/// @brief tf_ will be used internally by dockinghires. It will be initialized through the init_task_factory_ below
	core::pack::task::TaskFactoryOP tf_;
	/// @brief task_factory_ is used by outside movers to set the internal taskfactory. Does not actually override internal TF!
	/// init_task_factory_ instead acts as a starting point and the docking mover masks over init_task_factory
	core::Size repack_period_;

	bool norepack1_; // belongs to the packer task
	bool norepack2_; // belongs to the packer task

	// minimization
	core::Real min_tolerance_;
	bool nb_list_;
	std::string min_type_;
};

} // docking
} // protocols

#endif
