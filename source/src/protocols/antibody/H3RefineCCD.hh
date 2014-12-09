// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/H3RefineCCD.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody_H3RefineCCD_hh
#define INCLUDED_protocols_antibody_H3RefineCCD_hh


#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/loops/Loops.hh>

#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/ChangeFoldTreeMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/H3RefineCCD.fwd.hh>


using namespace core;
namespace protocols {
namespace antibody {

class H3RefineCCD: public moves::Mover {


public:

	/// @brief default constructor
	H3RefineCCD();

	/// @brief constructor with arguments
	H3RefineCCD(AntibodyInfoOP antibody_info, CDRNameEnum loop_name);

	/// @brief constructor with arguments
	H3RefineCCD(AntibodyInfoOP antibody_info, CDRNameEnum loop_name,
	            core::scoring::ScoreFunctionCOP highres_scorefxn );

	/// @brief constructor with arguments
	H3RefineCCD( AntibodyInfoOP antibody_info);


	void set_task_factory(pack::task::TaskFactoryCOP tf);

	virtual protocols::moves::MoverOP clone() const;

	/// @brief default destructor
	~H3RefineCCD();


	void pass_start_pose(core::pose::Pose & start_pose);
	void turn_on_benchmark() {
		benchmark_=true;
	}
	void turn_off_flank_relax() {
		flank_relax_ = false;
	}
	void turn_off_h3_filter() {
		H3_filter_=false;
	}

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void set_highres_score_func(core::scoring::ScoreFunctionCOP highres_scorefxn) {
		highres_scorefxn_ = highres_scorefxn->clone();
	}

private:

	AntibodyInfoOP ab_info_;
	CDRNameEnum loop_name_;
	loops::Loop the_loop_;
	core::Real high_cst_;
	core::Size loop_begin_,loop_end_,loop_size_;
	core::Size cutpoint_;
	core::Size n_small_moves_;
	core::Size inner_cycles_, outer_cycles_;
	bool include_neighbors_;

	core::kinematics::MoveMapOP cdrh3_map_, flank_cdrh3_map_;
	std::string minimization_type_;

	bool user_defined_;
	bool benchmark_;
	bool is_camelid_;


	// the objects
	moves::ChangeFoldTreeMoverOP change_FT_to_simpleloop_;
	moves::ChangeFoldTreeMoverOP change_FT_to_flankloop_;
	simple_moves::MinMoverOP loop_min_mover_;
	moves::SequenceMoverOP wiggle_cdr_h3_;
	simple_moves::PackRotamersMoverOP loop_repack_;
	moves::MonteCarloOP mc_;

	utility::vector1< bool> allow_repack_;

	core::Real init_temp_, last_temp_, gamma_;

	void set_default();
	void init();
	void finalize_setup( core::pose::Pose & pose );


	core::Real min_tolerance_;
	core::Real high_move_temp_;
	core::Real neighbor_dist_;

	// score functions
	core::scoring::ScoreFunctionOP highres_scorefxn_;


	/// @brief actually enables H3 filter for H3 operations
	bool H3_filter_;


	/// @brief number of flanking residues:default 5
	core::Size flank_size_;

	/// @brief relax flanking regions of h3
	bool flank_relax_;

	core::pose::Pose start_pose_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;

	/// @brief just refine input loop
	bool refine_input_loop_;

	core::Size max_cycle_close_trial_;

};







} // namespace antibody
} // namespace protocols

#endif
