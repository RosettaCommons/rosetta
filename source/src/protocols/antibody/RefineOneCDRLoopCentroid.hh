// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/RefineOneCDRLoopCentroid.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_RefineOneCDRLoopCentroid_hh
#define INCLUDED_protocols_antibody_RefineOneCDRLoopCentroid_hh


#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/RefineOneCDRLoopCentroid.fwd.hh>


namespace protocols {
namespace antibody {

class RefineOneCDRLoopCentroid: public moves::Mover {


public:

	/// @brief constructor with arguments
	RefineOneCDRLoopCentroid(AntibodyInfoCOP antibody_info,
		CDRNameEnum const & loop_name);

	RefineOneCDRLoopCentroid(AntibodyInfoCOP antibody_info,
		CDRNameEnum const & loop_name,
		core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief constructor with arguments
	RefineOneCDRLoopCentroid( loops::Loop const & a_cdr_loop);

	RefineOneCDRLoopCentroid( loops::Loop const & a_cdr_loop,
		core::scoring::ScoreFunctionCOP scorefxn );

	/// @brief default destructor
	~RefineOneCDRLoopCentroid();


	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;


	void set_benchmark(bool const & setting) {
		benchmark_ = setting;
	}
	void set_snugfit(bool const & setting) {
		snug_fit_ = setting;
	}
	void set_refine_input_loop(bool const & setting) {
		refine_input_loop_ = setting;
	}

	void set_score_function(core::scoring::ScoreFunctionCOP scorefxn);


private:
	void set_default();
	void finalize_setup( core::pose::Pose const & pose );
	void loop_centroid_relax(
		core::pose::Pose & pose,
		Size const loop_begin,
		Size const loop_end );

private:
	loops::Loop the_cdr_loop_;

	bool benchmark_;
	bool snug_fit_;
	bool refine_input_loop_;
	core::scoring::ScoreFunctionOP lowres_scorefxn_;

};


} // namespace antibody
} // namespace protocols

#endif
