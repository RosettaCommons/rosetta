// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file docking_initialization_protocols
/// @brief initialization protocols for docking
/// @detailed
/// @author  Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockingHiRes_hh
#define INCLUDED_protocols_symmetric_docking_SymDockingHiRes_hh

// Package headers
#include <protocols/symmetric_docking/SymDockingHiRes.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.fwd.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.fwd.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.fwd.hh>

// For symmetry
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

class SymDockingHiRes : public moves::Mover
{
	typedef core::Real Real;
	typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;


public:

	// default constructor
	SymDockingHiRes();

	/* // constructor with arguments
	SymDockingHiRes(
		core::scoring::ScoreFunctionOP scorefxn_in
	); */

	// constructor with arguments
	SymDockingHiRes(
		core::scoring::ScoreFunctionOP scorefxn_in,
		core::scoring::ScoreFunctionOP scorefxn_pack_in
	);

	moves::MoverOP clone() const;

	virtual ~SymDockingHiRes();

	void set_default( core::pose::Pose & pose );
	void set_move_map(core::kinematics::MoveMapOP movemap_in);
	void set_min_type( std::string min_type_in );
	void set_repack( bool repack_switch);
	void set_protocol( core::pose::Pose & pose );

	void set_dock_min_protocol();
	void set_dock_mcm_protocol( core::pose::Pose & pose );
	void set_dock_ppk_protocol( core::pose::Pose & pose );

	//void define_loops( core::pose::Pose const & pose, loops::Loops & loop_set, Real & interface_dist );

	moves::MonteCarloOP get_mc();

	// @brief turns on design of partner2 during docking. Experimental!
	void design( bool const des );
	bool design() const;

	// protocol functions
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	//void call_pack();  // Undefined function, commenting out to make python bindginds

	//void dock_mcm_protocol( core::pose::Pose & pose );

	/*void classic_mcm_protocol(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scorefxn,
		protocols::moves::MonteCarloOP monteCarlo,
		core::Size num_cycles,
		core::Size repack_every_Nth
	); */

	void task_factory( core::pack::task::TaskFactoryOP task );
	core::pack::task::TaskFactoryOP & task_factory();

private:
	// protocol stuff
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionCOP scorefxn_pack_;
	core::kinematics::MoveMapOP movemap_;
	simple_moves::symmetry::SymMinMoverOP min_mover_;
	moves::SequenceMoverOP docking_highres_protocol_mover_;
	moves::MonteCarloOP mc_;
	simple_moves::symmetry::SymPackRotamersMoverOP pack_interface_repack_;
	simple_moves::symmetry::SymRotamerTrialsMoverOP pack_rottrial_;

	// docking
	core::Real trans_magnitude_, rot_magnitude_;
	bool chi_, bb_;
	bool repack_switch_; // master switch to turn repacking during docking on/off. Only makes sense if repack_period > 0
	bool design_;
	bool rtmin_, scmin_;

	// packing
	/// @brief utility function to set up packer options for internal task factory
	void setup_packing( core::pose::Pose & pose );
	/// @brief tf_ will be used internally by dockinghires. It will be initialized through the init_task_factory_ below
	core::pack::task::TaskFactoryOP tf_;
	/// @brief task_factory_ is used by outside movers to set the internal taskfactory. Does not actually override internal TF!
	/// init_task_factory_ instead acts as a starting point and the docking mover masks over init_task_factory
	core::pack::task::TaskFactoryOP init_task_factory_;
	core::Size repack_period_;
	core::Real temperature_;

	// minimization
	core::Real min_tolerance_;
	bool nb_list_;
	std::string min_type_;

};

} // symmetrical_docking
} // protocols

#endif
