// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/RefineOneCDRLoop.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_RefineOneCDRLoop_hh
#define INCLUDED_protocols_antibody_RefineOneCDRLoop_hh


#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/RefineOneCDRLoop.fwd.hh>


namespace protocols {
namespace antibody {

class RefineOneCDRLoop: public moves::Mover {


public:

	/// @brief default constructor
	RefineOneCDRLoop();

	/// @brief constructor with arguments
	RefineOneCDRLoop( AntibodyInfoOP antibody_info);

	/// @brief constructor with arguments
	RefineOneCDRLoop(AntibodyInfoOP antibody_info, std::string refine_mode);

	/// @brief constructor with arguments
	RefineOneCDRLoop(AntibodyInfoOP antibody_info, std::string refine_mode,
		core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief constructor with arguments
	RefineOneCDRLoop(AntibodyInfoOP antibody_info,
		CDRNameEnum const & cdr_loop_name,
		std::string refine_mode,
		core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief default destructor
	~RefineOneCDRLoop();


	void turn_on_benchmark() {
		benchmark_=true;
	}
	void set_h3_filter(bool setting) {
		H3_filter_=setting;
	}
	void set_num_filter_tries(core::Size setting) {
		num_filter_tries_=setting;
	}
	void set_flank_relax(bool setting) {
		flank_relax_=setting;
	}
	void set_flank_size(core::Size setting) {
		flank_size_=setting;
	}
	void set_refine_mode(std::string refine_mode) {
		refine_mode_ = refine_mode;
	}
	virtual void apply( core::pose::Pose & pose );


	void set_score_function(core::scoring::ScoreFunctionCOP scorefxn);

	void pass_start_pose(core::pose::Pose & start_pose);
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;


private:

	AntibodyInfoOP ab_info_;
	bool user_defined_;
	bool benchmark_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Size flank_size_;
	bool flank_relax_;
	core::pose::Pose start_pose_;
	core::pack::task::TaskFactoryOP tf_;
	bool H3_filter_;
	core::Size num_filter_tries_;
	std::string refine_mode_;
	CDRNameEnum cdr_loop_name_;
	core::Real high_cst_;

	void set_default();
	void init();
	void finalize_setup( core::pose::Pose & pose );

};


} // namespace antibody
} // namespace protocols

#endif


