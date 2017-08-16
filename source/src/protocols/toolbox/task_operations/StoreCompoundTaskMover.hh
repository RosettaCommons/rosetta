// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/StoreCompoundTaskMover.cc
/// @brief  Combine tasks using boolean logic for residues that are packable or designable, as
/// well as for residue specific AA sets, and store the resulting task in the pose's cacheable data.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_StoreCompoundTaskMover_hh
#define INCLUDED_protocols_toolbox_task_operations_StoreCompoundTaskMover_hh

//unit headers
#include <protocols/toolbox/task_operations/StoreCompoundTaskMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class StoreCompoundTaskMover : public protocols::moves::Mover {

public:

	typedef std::vector< std::pair< core::pack::task::PackerTaskOP, boolean_operations > > CompoundTask;
	typedef CompoundTask::iterator task_iterator;
	typedef CompoundTask::const_iterator const_task_iterator;

	typedef std::vector< std::pair< core::pack::task::TaskFactoryOP, boolean_operations > > CompoundFactory;
	typedef CompoundFactory::iterator factory_iterator;
	typedef CompoundFactory::const_iterator const_factory_iterator;

	StoreCompoundTaskMover();
	~StoreCompoundTaskMover();

	// Boolean logic functions for the different modes.
	void CompoundPackableTask( core::pose::Pose const & pose, core::Size & total_residue, core::pack::task::PackerTaskOP & task );
	void CompoundDesignableTask( core::pose::Pose const & pose, core::Size & total_residue, core::pack::task::PackerTaskOP & task );
	// void CompoundAminoAcidSetTask( core::Size & total_residue, core::pack::task::PackerTaskOP & task ); TODO (Jacob)

	void apply( core::pose::Pose & pose  ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void task_clear();
	task_iterator task_begin();
	const_task_iterator task_begin() const;
	task_iterator task_end();
	const_task_iterator task_end() const;
	void factory_clear();
	factory_iterator factory_begin();
	const_factory_iterator factory_begin() const;
	factory_iterator factory_end();
	const_factory_iterator factory_end() const;
	void task_name( std::string const & tn );
	void mode( std::string const & md );
	void true_behavior( std::string const & tb );
	void false_behavior( std::string const & fb );
	void invert( bool const inv );
	void verbose( bool const verb );
	void overwrite( bool const ow );

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose
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
	CompoundTask compound_task_;
	CompoundFactory compound_factory_;
	std::string task_name_;
	std::string mode_;
	std::string true_behavior_;
	std::string false_behavior_;
	bool invert_;
	bool verbose_;
	bool overwrite_;
};


} // task_operations
} // toolbox
} // protocols

#endif
