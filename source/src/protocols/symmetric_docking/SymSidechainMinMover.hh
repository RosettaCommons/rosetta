// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sid Chaudhury ( carbon copy by Ingemar André )
#ifndef INCLUDED_protocols_symmetric_docking_SymSidechainMinMover_hh
#define INCLUDED_protocols_symmetric_docking_SymSidechainMinMover_hh

// Unit headers
#include <protocols/symmetric_docking/SymSidechainMinMover.fwd.hh>

// Project headers
#include <protocols/scoring/Interface.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

class SymSidechainMinMover : public moves::Mover {
public:

	// default constructor
	SymSidechainMinMover();

	//constructors with arguments
	SymSidechainMinMover( core::scoring::ScoreFunctionCOP scorefxn_in );
	SymSidechainMinMover( core::scoring::ScoreFunctionCOP scorefxn_in,core::kinematics::MoveMapOP movemap_in );
	SymSidechainMinMover( core::scoring::ScoreFunctionCOP scorefxn_in, core::pack::task::PackerTaskOP);
	SymSidechainMinMover( core::scoring::ScoreFunctionCOP scorefxn_in, core::pack::task::TaskFactoryOP );

	// destructor
	~SymSidechainMinMover();

	void set_default_options();

	void set_minmover( protocols::simple_moves::MinMoverOP minmover_in );

	void update_movemap( core::pose::Pose & pose );

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

protected:

	core::scoring::ScoreFunctionCOP scorefxn_;
	protocols::simple_moves::MinMoverOP minmover_;
	core::kinematics::MoveMapOP movemap_;
	core::pack::task::PackerTaskOP task_;
	core::pack::task::TaskFactoryOP tf_;
	bool update_movemap_;

};


class SymInterfaceSidechainMinMover : public SymSidechainMinMover {
public:

	//default constructor
	SymInterfaceSidechainMinMover();

	// constructor with arguments
	SymInterfaceSidechainMinMover(
		core::scoring::ScoreFunctionCOP scorefxn_in,
		core::Real interface_dist_in=8.0
	);

	// destructor
	~SymInterfaceSidechainMinMover();

	void set_default_options();

	void set_interface_dist( core::Real interface_dist_in);

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	core::Real interface_dist_;

	protocols::scoring::InterfaceOP interface_;

};

} //symmetric_docking
} // protocols


#endif
