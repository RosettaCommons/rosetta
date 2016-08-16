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


#ifndef INCLUDED_protocols_docking_DockMCMProtocol_hh
#define INCLUDED_protocols_docking_DockMCMProtocol_hh

// Unit Headers
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/docking/DockFilters.fwd.hh>

#include <protocols/docking/DockMCMCycle.fwd.hh>
#include <utility/vector1.hh>
#include <core/kinematics/MoveMap.fwd.hh>

namespace protocols {
namespace docking {

class DockMCMProtocol : public DockingHighRes
{
	typedef core::Real Real;
public:

	/// @brief Default constructor
	DockMCMProtocol();

	/// @brief Constructor with one argument - the jump number.
	DockMCMProtocol(
		core::Size const rb_jump
	);

	/// @brief Constructor with three arguments. The first is the jump number, the second is the docking scorefunction
	/// and the third is the packing scorefxn.
	DockMCMProtocol(
		core::Size const rb_jump,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	/// @brief Constructor with two arguments. The first is the DockJumps, the second is a scorefunction that will be
	/// used for docking and packing.
	DockMCMProtocol(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn
	);

	/// @brief Constructor with three arguments. The first is the DockJumps, the second is the docking scorefunction
	/// and the third is the packing scorefunction.
	DockMCMProtocol(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	// destructor
	~DockMCMProtocol();

	//clone
	protocols::moves::MoverOP clone() const;

	/// @brief Performs the setup specific to this subclass of DockingHighRes (sets up filters).
	void init();

	void set_filter( DockingHighResFilterOP filter );

	// protocol functions
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void set_move_map(core::kinematics::MoveMapOP movemap );
	void set_second_cycle(Size const & num);
	void set_first_cycle(Size const & num);

	core::scoring::ScoreFunctionCOP scorefxn_docking() const;
	core::scoring::ScoreFunctionCOP scorefxn_packing() const;
	friend std::ostream & operator<<(std::ostream& out, const DockMCMProtocol & dmp );

private:
	DockingHighResFilterOP filter_;
	DockMCMCycleOP dock_mcm_; //JQX: make it as a memmber
	bool movemap_reset_ ;
	core::kinematics::MoveMapOP movemap_;
	Size num_of_first_cycle_;
	Size num_of_second_cycle_;
};

} // docking
} // protocols

#endif
