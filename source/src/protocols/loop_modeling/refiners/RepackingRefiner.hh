// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_refiners_RepackingRefiner_HH
#define INCLUDED_protocols_loop_modeling_refiners_RepackingRefiner_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

/// @brief Refine sampled loops using sidechain repacking.
class RepackingRefiner : public LoopMover {

public:

	/// @brief Default constructor.
	RepackingRefiner(Size repack_period=1);

	/// @copydoc LoopMover::get_name
	string get_name() const { return "RepackingRefiner"; }

	/// @copydoc LoopMover::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose);

	/// @brief Get the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP score_function);

	/// @brief Get the task factory to be used on the next call to apply().
	/// @details If no task factory has been set, this will raise an exception.
	core::pack::task::TaskFactoryOP get_task_factory();

	/// @brief Get the task factory to be used on the next call to apply().
	/// @details If no task factory has been set, the fallback will be returned.
	core::pack::task::TaskFactoryOP get_task_factory(
		core::pack::task::TaskFactoryOP fallback);

	/// @brief Set the task factory to be used on the next call to apply().
	void set_task_factory(core::pack::task::TaskFactoryOP task_factory);

	/// @brief Return how often this refiner actually repacks the pose.
	Size get_repack_period() const;

	/// @brief Specify how often this refiner should actually repack the pose.
	void set_repack_period(Size period);

protected:

	/// @brief Repack the pose within 10A of the loops being sampled.
	bool do_apply(Pose & pose);

private:

	Size repack_period_;
	Size iteration_counter_;

};

}
}
}

#endif
