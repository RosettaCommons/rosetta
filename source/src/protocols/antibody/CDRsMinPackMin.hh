// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/CDRsMinPackMin.hh
/// @brief Main CDR minimization method in protocols/antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody_CDRsMinPackMin_hh
#define INCLUDED_protocols_antibody_CDRsMinPackMin_hh

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

#include <protocols/antibody/ModelCDRH3.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/CDRsMinPackMin.fwd.hh>
#include <protocols/antibody/RefineBetaBarrel.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace antibody {

class CDRsMinPackMin: public moves::Mover {
public:

	CDRsMinPackMin();

	CDRsMinPackMin(AntibodyInfoOP antibody_info);

	CDRsMinPackMin(
		AntibodyInfoOP antibody_info,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf,
		core::kinematics::MoveMapOP movemap
	);

	// default destructor
	~CDRsMinPackMin() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void set_task_factory(core::pack::task::TaskFactoryCOP tf);
	void set_move_map(core::kinematics::MoveMapCOP movemap);

	// simple inline setters
	void set_sc_min (bool scmin) {
		sc_min_ = scmin ;
	}
	void set_rt_min (bool rtmin) {
		rt_min_ = rtmin ;
	}
	void set_turnoff_minimization(bool setting) {
		turnoff_minimization_=setting;
	}
	void set_alter_foldtree(bool setting) {
		alter_foldtree_ = setting;
	}

	void show( std::ostream & out=std::cout ) const override;
	friend std::ostream & operator<<(std::ostream& out, const CDRsMinPackMin & ab_m_2 );


private:
	bool sc_min_;
	bool rt_min_;
	bool turnoff_minimization_;

	// do not touch the foldtree
	bool alter_foldtree_;

	// Benchmark mode for shorter_cycles
	bool benchmark_;
	core::Size update_rounds_;
	//to update the task factory and movemap auto-ly

	void finalize_setup( core::pose::Pose & pose );
	void init();
	bool user_defined_; // for constructor options passed to init

	AntibodyInfoOP ab_info_;
	core::scoring::ScoreFunctionOP loop_scorefxn_highres_;
	protocols::moves::SequenceMoverOP cdr_sequence_move_ ;
	core::kinematics::MoveMapOP allcdr_map_;
	core::pack::task::TaskFactoryOP tf_;
	std::string min_type_;
	core::Real Temperature_;
	core::Real min_tolerance_;


}; // class CDRsMinPackMin


} // namespace antibody
} // namespace protocols

#endif

