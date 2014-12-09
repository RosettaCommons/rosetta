// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/LHRepulsiveRamp.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody_LHRepulsiveRamp_hh
#define INCLUDED_protocols_antibody_LHRepulsiveRamp_hh


#include <protocols/antibody/LHRepulsiveRamp.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/docking/types.hh>


using namespace core;
namespace protocols {
namespace antibody {

class LHRepulsiveRamp: public moves::Mover {


public:

	/// @brief default constructor
	LHRepulsiveRamp();

	/// @brief constructor with arguments

	LHRepulsiveRamp( docking::DockJumps const movable_jumps,
	                 core::scoring::ScoreFunctionCOP dock_scorefxn,
	                 core::scoring::ScoreFunctionCOP pack_scorefxn );

	virtual protocols::moves::MoverOP clone() const;

	/// @brief default destructor
	~LHRepulsiveRamp();

	void set_default();

	void set_dock_score_func(scoring::ScoreFunctionCOP dock_scorefxn ) {
		dock_scorefxn_ = dock_scorefxn->clone();
	}

	void set_pack_score_func(scoring::ScoreFunctionCOP pack_scorefxn) {
		pack_scorefxn_ = pack_scorefxn->clone();
	}

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;


	void set_task_factory(pack::task::TaskFactoryCOP tf);
	void set_move_map(kinematics::MoveMapCOP movemap);
	void set_dock_jump(docking::DockJumps jump);
	Real set_rot_mag  (core::Real rot_mag)  {
		return rot_mag_  =rot_mag;
	}
	Real set_trans_mag(core::Real trans_mag) {
		return trans_mag_=trans_mag;
	}


	void set_sc_min(bool sc_min) {
		sc_min_ = sc_min;
	}

	void set_rt_min(bool rt_min) {
		rt_min_ = rt_min;
	}

private:


	bool user_defined_;
	bool benchmark_;
	core::Size rep_ramp_cycles_;
	core::Real rot_mag_;
	core::Real trans_mag_;
	core::Size num_repeats_;

	scoring::ScoreFunctionOP dock_scorefxn_;
	scoring::ScoreFunctionOP pack_scorefxn_;

	void init();

	void repulsive_ramp( pose::Pose & pose_in, loops::Loops loops_in );


	//packer task
	docking::DockJumps jump_;
	pack::task::TaskFactoryOP tf_;
	kinematics::MoveMapOP movemap_;
	bool sc_min_;
	bool rt_min_;
};







} // namespace antibody
} // namespace protocols

#endif








