// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sid Chaudhury
#ifndef INCLUDED_protocols_docking_SidechainMinMover_hh
#define INCLUDED_protocols_docking_SidechainMinMover_hh

// Unit headers
#include <protocols/docking/SidechainMinMover.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>

// Project headers
#include <protocols/scoring/Interface.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/PackerTask.hh>

#include <protocols/simple_moves/MinMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace docking {

class SidechainMinMover : public DockingHighRes {
public:

	/// @brief Default constructor
	SidechainMinMover();

	/// @brief Constructor with one argument - the scoerfunction to minimize with
	SidechainMinMover( core::scoring::ScoreFunctionOP scorefxn );

	/// @brief Constructor with two arguments.  The first is the scorefunction to minimize with, the second is a movemap
	SidechainMinMover( core::scoring::ScoreFunctionOP scorefxn,core::kinematics::MoveMapOP movemap );

	/// @brief Constructor with two arguments. The first is the scorefunction to minimize with, the second is a task
	SidechainMinMover( core::scoring::ScoreFunctionOP scorefxn, core::pack::task::PackerTaskOP);

	/// @brief Constructor with two arguments. The first is the scorefunction to minimize with, the second is a taskfactory
	SidechainMinMover( core::scoring::ScoreFunctionOP scorefxn, core::pack::task::TaskFactoryCOP );

	/// @brief Constructor with two arguments.  The first is the jump that docking will occur over, the second is the
	///  scorefunction to minimize with.
	SidechainMinMover( core::Size rb_jump, core::scoring::ScoreFunctionOP scorefxn );

	// destructor
	~SidechainMinMover();


	/// @brief Completes the setup of a default instantiation of a SideChainMinMover
	void set_default();


	void set_minmover( protocols::simple_moves::MinMoverOP minmover );

	void update_movemap( core::pose::Pose & pose );

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

protected:
	protocols::simple_moves::MinMoverOP minmover_;
	core::kinematics::MoveMapOP movemap_;
	core::pack::task::PackerTaskOP task_;
	bool update_movemap_;
};


class InterfaceSidechainMinMover : public SidechainMinMover {
public:

	//default constructor
	InterfaceSidechainMinMover();

	// constructor with arguments
	InterfaceSidechainMinMover(
		core::Size rb_jump,
		core::scoring::ScoreFunctionOP scorefxn,
		core::Real interface_dist=8.0
	);

	// destructor
	~InterfaceSidechainMinMover();

	void set_default();

	void set_interface_dist( core::Real interface_dist);
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Real interface_dist_;

	protocols::scoring::InterfaceOP interface_;
};

} //docking
} // protocols


#endif
