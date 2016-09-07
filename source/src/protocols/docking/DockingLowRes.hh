// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockinLowRes
/// @brief low resolution mode for docking
/// @details
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov


#ifndef INCLUDED_protocols_docking_DockingLowRes_hh
#define INCLUDED_protocols_docking_DockingLowRes_hh

#include <protocols/docking/types.hh>
#include <protocols/docking/DockingLowRes.fwd.hh>


// Package headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

#include <string>

#include <utility/vector1.hh>
#include <iostream>

#ifdef PYROSETTA
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#endif


namespace protocols {
namespace docking {

class DockingLowRes : public moves::Mover
{
	//typedef core::Real Real;
public:

	/// @brief Default constructor
	DockingLowRes();

	// destructor
	~DockingLowRes() override;

	/// @brief Constructor with two arguments.  The first is scorefunction to be used for docking, the second is the
	///  jump to dock over.
	DockingLowRes(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size const rb_jump=1
	);

	/// @brief Constructor with two arguments.  The first is scorefunction to be used for docking, the second is the
	///  DockJumps.
	DockingLowRes(
		core::scoring::ScoreFunctionCOP scorefxn,
		DockJumps const movable_jumps
	);

	protocols::moves::MoverOP clone() const override;

	/// @brief Assigns default values to primitive members
	void set_default();
	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	moves::MonteCarloOP get_mc();

	// protocol functions
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	/// @brief Perform a cycle of rigid-body Monte Carlo moves
	void rigid_body_trial( core::pose::Pose & pose );

	// option setters
	void set_inner_cycles( core::Size inner_cycles ) { inner_cycles_=inner_cycles; }
	void set_outer_cycles( core::Size outer_cycles ) { outer_cycles_=outer_cycles; }

	void set_trans_magnitude( core::Real trans_magnitude);
	void set_rot_magnitude( core::Real rot_magnitude);

	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, const DockingLowRes & dp );

	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	// Add by dK
	// docking
	DockJumps movable_jumps_;

	// Add by DK
	protocols::moves::SequenceMoverOP docking_lowres_protocol_;
	core::scoring::ScoreFunctionCOP scorefxn_;

protected:
	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	virtual void finalize_setup( core::pose::Pose & pose);

private:
	// protocol stuff
	//core::scoring::ScoreFunctionCOP scorefxn_;
	core::kinematics::MoveMapOP movemap_;
	protocols::rigid::RigidBodyPerturbNoCenterMoverOP rb_mover_;
	protocols::moves::MonteCarloOP mc_;

	// Comment out by DK
	// docking
	//DockJumps movable_jumps_;
	core::Real trans_magnitude_, rot_magnitude_, accept_rate_;
	core::Size inner_cycles_, outer_cycles_; //rb_jump_
	bool chi_, bb_, nb_list_;
	core::Real temperature_;

	/// @brief Sets up the instance of DockingLowRes and initializes all members based on values passed in at construction
	///  or via the command line.
	void init(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionCOP scorefxn
	);
};

} // docking
} // protocols

#endif
