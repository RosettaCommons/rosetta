// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/RefineBetaBarrel.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_RefineBetaBarrel_hh
#define INCLUDED_protocols_antibody_RefineBetaBarrel_hh

#include <protocols/moves/Mover.hh>
#include <protocols/antibody/RefineBetaBarrel.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody/LHRepulsiveRamp.fwd.hh>
#include <protocols/antibody_legacy/LHSnugFitLegacy.fwd.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>
#include <protocols/docking/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


namespace protocols {
namespace antibody {


class RefineBetaBarrel: public moves::Mover {


public:
	/// @brief default constructor
	RefineBetaBarrel();

	/// @brief default destructor
	~RefineBetaBarrel();

	RefineBetaBarrel(AntibodyInfoOP antibody_info);

	RefineBetaBarrel(AntibodyInfoOP antibody_info,
		core::scoring::ScoreFunctionCOP dock_scorefxn,
		core::scoring::ScoreFunctionCOP pack_scorefxn);

	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	void set_task_factory(core::pack::task::TaskFactoryCOP tf);


	void set_dock_score_func(core::scoring::ScoreFunctionCOP dock_scorefxn );

	void set_pack_score_func(core::scoring::ScoreFunctionCOP pack_scorefxn);

	void turn_off_repulsive_ramp() {
		repulsive_ramp_ = false;
	}


	void set_sc_min(bool sc_min) {
		sc_min_ = sc_min;
	}

	void set_rt_min(bool rt_min) {
		rt_min_ = rt_min;
	}

private:
	bool sc_min_;
	bool rt_min_;
	bool user_defined_;
	bool repulsive_ramp_;
	AntibodyInfoOP ab_info_;
	core::pack::task::TaskFactoryOP tf_;


	void init( );
	void finalize_setup(core::pose::Pose & pose_in );

	LHRepulsiveRampOP lh_repulsive_ramp_;
	LHSnugFitLegacyOP lh_snugfit_;
	docking::DockMCMProtocolOP dock_mcm_protocol_;

	core::scoring::ScoreFunctionOP dock_scorefxn_;
	core::scoring::ScoreFunctionOP pack_scorefxn_;

	core::kinematics::MoveMapOP cdr_dock_map_;

	docking::DockJumps LH_dock_jump_;
};


}// namespace antibody
}// namespace protocols


#endif


