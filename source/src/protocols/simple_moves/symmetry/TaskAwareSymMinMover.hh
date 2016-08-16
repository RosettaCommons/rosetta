// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/symmetry/TaskAwareSymMinMover.hh
/// @brief  Minimize only those residues specified by a set of TaskOperations in a symmetric pose
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_symmetry_TaskAwareSymMinMover_hh
#define INCLUDED_protocols_simple_moves_symmetry_TaskAwareSymMinMover_hh

// Unit Headers
#include <protocols/simple_moves/symmetry/TaskAwareSymMinMover.fwd.hh>
#include <protocols/simple_moves/symmetry/TaskAwareSymMinMoverCreator.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/types.hh>

// Utility Headers

// C++ Headers

namespace protocols {
namespace simple_moves {
namespace symmetry {

//class TaskAwareSymMinMover : public protocols::moves::Mover {
class TaskAwareSymMinMover : public protocols::simple_moves::TaskAwareMinMover { //Jeliazkov experiment
public:

	// default constructor
	TaskAwareSymMinMover();

	// constructor with arguments
	TaskAwareSymMinMover(
		protocols::simple_moves::MinMoverOP minmover_in,
		core::pack::task::TaskFactoryCOP factory_in
	);


	TaskAwareSymMinMover(const TaskAwareSymMinMover& rval);

	virtual std::string get_name() const { return "TaskAwareSymMinMover"; }
	virtual void apply( core::pose::Pose & pose );

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose );

private:

	std::string scorefxn_name_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool min_chi_;
	bool min_bb_;
	bool min_rb_; //jump
	std::string min_type_;
	core::Real tolerance_;
	//core::pack::task::TaskFactoryOP task_factory_;
	/// @brief OP for MinMover
	protocols::simple_moves::MinMoverOP minmover_;
	/// @brief OP for constant task factory for nonconstant tasks, if present
	core::pack::task::TaskFactoryCOP factory_;
	bool designable_only_;

};

} //namespace symmetry
} //namespace simple_moves
} //namespace protocols

#endif
