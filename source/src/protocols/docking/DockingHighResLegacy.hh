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


#ifndef INCLUDED_protocols_docking_DockingHighResLegacy_hh
#define INCLUDED_protocols_docking_DockingHighResLegacy_hh

//Unit headers
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockingHighResLegacy.fwd.hh>

// Package headers

//Project headers
#include <protocols/loops/Loops.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace docking {

/// @brief this mover does the high resolution refinement stage of the RosettaDock algorithm
class DockingHighResLegacy : public DockingHighRes
{
	typedef core::Real Real;
public:

	// default constructor
	DockingHighResLegacy();

	// constructor with arguments
	DockingHighResLegacy(
		int rb_jump,
		core::scoring::ScoreFunctionOP scorefxn

	);

	// constructor with arguments
	DockingHighResLegacy(
		int rb_jump,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	// constructor with arguments
	DockingHighResLegacy(
		DockJumps const movable_jumps,
		core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreFunctionOP scorefxn_pack
	);

	~DockingHighResLegacy() override;

	//clone
	protocols::moves::MoverOP clone() const override;

	void set_default( core::pose::Pose & pose );
	void set_move_map(core::kinematics::MoveMapOP movemap_in);
	void set_min_type( std::string min_type_in );
	void set_repack( bool repack_switch);
	void set_protocol( core::pose::Pose & pose );
	virtual void set_task_factory( core::pack::task::TaskFactoryOP task );

	void define_loops( core::pose::Pose const & pose, loops::LoopsOP loop_set, Real & interface_dist );

	void set_trans_magnitude( core::Real trans_magnitude);
	void set_rot_magnitude( core::Real rot_magnitude);
	void set_rt_min( bool rt_min ){ rt_min_ = rt_min; }
	bool rt_min() { return rt_min_; }
	void set_sc_min( bool sc_min ){ sc_min_ = sc_min; }
	bool sc_min() { return sc_min_; }

	moves::MonteCarloOP get_mc();
	core::Real get_trans_magnitude() {return trans_magnitude_;}
	core::Real get_rot_magnitude() {return rot_magnitude_;}

	// protocol functions
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void set_dock_min_protocol();
	void set_dock_mcm_protocol( core::pose::Pose & pose );
	void set_dock_ppk_protocol( core::pose::Pose & pose );

	// @brief turns on design of partner2 during docking. Experimental!
	void design( bool const des );
	bool design() const;

private:
	core::kinematics::MoveMapOP movemap_;
	moves::SequenceMoverOP docking_highres_protocol_mover_;
	moves::MonteCarloOP mc_;

	// docking
	core::Real trans_magnitude_, rot_magnitude_;
	bool chi_, bb_;
	bool repack_switch_; // master switch to turn repacking during docking on/off. Only makes sense if repack_period > 0
	bool design_;
	bool rt_min_;
	bool sc_min_;

	// packing
	/// @brief utility function to set up packer options for internal task factory
	void setup_packing( core::pose::Pose & pose );
	core::pack::task::TaskFactoryOP init_task_factory_;

	core::Size repack_period_;
	core::Real temperature_;

	// minimization
	core::Real min_tolerance_;
	bool nb_list_;
	std::string min_type_;
};

} // docking
} // protocols

#endif
