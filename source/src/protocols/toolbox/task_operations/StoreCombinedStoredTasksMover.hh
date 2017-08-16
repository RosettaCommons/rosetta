// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/StoreCombinedStoredTasksMover.hh
/// @brief Headers for StoreCombinedStoredTasksMover class -- DEPRECATED (see .cc file)
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_StoreCombinedStoredTasksMover_hh
#define INCLUDED_protocols_toolbox_task_operations_StoreCombinedStoredTasksMover_hh

//unit headers
#include <protocols/toolbox/task_operations/StoreCombinedStoredTasksMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @brief mover that can be used to save or restore a task at an arbitrary
/// point during a rosetta scripts protocol. other task operations, movers,
/// or filters can be set up to access tasks saved by this mover during their
/// apply calls.
class StoreCombinedStoredTasksMover : public protocols::moves::Mover {

public:

	StoreCombinedStoredTasksMover();
	~StoreCombinedStoredTasksMover();

	void apply( core::pose::Pose & pose  ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string task1_, task2_, task_name_, operator_;
	bool overwrite_, invert_;
};


} // task_operations
} // toolbox
} // protocols

#endif

