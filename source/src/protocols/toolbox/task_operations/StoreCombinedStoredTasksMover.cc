// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/StoreCombinedStoredTasksMover.cc
/// @brief The StoreCombinedStoredTasksMover HAS BEEN DEPRECATED. It was used to combine
/// several STMStoredTasks (protocols/toolbox/task_operations/STMStoredTask.cc), but has been
/// replaced by the StoreCompoundTaskMover
/// (protocols/toolbox/task_operations/StoreCompoundTaskMover.cc), which is much more flexible.
/// It is being moved out of devel/matdes/ into protocols/toolbox/task_operations/ to support
/// the publication King et al., Nature 2014.
/// @author Jacob Bale (balej@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/StoreCombinedStoredTasksMover.hh>
#include <protocols/toolbox/task_operations/StoreCombinedStoredTasksMoverCreator.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <protocols/toolbox/task_operations/STMStoredTask.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.StoreCombinedStoredTasksMover" );

namespace protocols {
namespace toolbox {
namespace task_operations {

// @brief default constructor
StoreCombinedStoredTasksMover::StoreCombinedStoredTasksMover() {}

// @brief destructor
StoreCombinedStoredTasksMover::~StoreCombinedStoredTasksMover() {}

void
StoreCombinedStoredTasksMover::apply( core::pose::Pose & pose )
{
	// Create an initial task to subsequently be restricted based on the two preivously stored tasks.
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	// Retrieve the user-defined stored tasks from the pose.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		utility_exit_with_message("Your pose does not have CacheableData of type STM_STORED_TASKS");
	} else {
		core::Size total_residue;
		// Grab a reference to the data and check that the user-defined tasks are present.
		protocols::toolbox::task_operations::STMStoredTask & stored_tasks = *( utility::pointer::static_pointer_cast< protocols::toolbox::task_operations::STMStoredTask > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) );
		if ( !stored_tasks.has_task(task1_) ) {
			utility_exit_with_message("No stored task with the name " + task1_ + " found");
		} else if ( !stored_tasks.has_task(task2_) ) {
			utility_exit_with_message("No stored task with the name " + task2_ + " found");
		} else {
			// Only need to go through the asymmetric unit if the pose is symmetric.
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
				total_residue = symm_info->num_independent_residues();
			} else {
				total_residue = pose.total_residue();
			}
		}
		// Loop throught the residues and determine the state of each residue for each retrieved task.
		// The logic applies to the packable positions.
		// Ex. if invert=false and operator=AND, then only those residues that were packable in task1 AND task2 will be packable in the new task.
		// Ex. if invert=true and operator=AND, then only those residues that were not packable in task1 AND task2 will be packable in the new task.
		std::string select_packable_pos("select design_positions, resi ");
		for ( core::Size resi=1; resi<=total_residue; ++resi ) {
			// If taking the intersection, then prevent repacking at all residues except those that were packable in both task1 AND task2.
			if ( !invert_ && operator_ == "AND" ) {
				if ( !(stored_tasks.get_task( task1_ )->being_packed( resi ) && stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
				// If taking the inverse of the intersection, then prevent repacking at all residues that were packable in both task1 AND task2.
			} else if ( invert_ && operator_ == "AND" ) {
				if ( (stored_tasks.get_task( task1_ )->being_packed( resi ) && stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
				// If taking the union, then prevent repacking at all residues except those that were packable in task1 OR task2.
			} else if ( !invert_ && operator_ == "OR" ) {
				if ( !(stored_tasks.get_task( task1_ )->being_packed( resi ) || stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
				// If taking the inverse of the union, then prevent repacking at all residues that were packable in task1 OR task2.
			} else if ( invert_ && operator_ == "OR" ) {
				if ( (stored_tasks.get_task( task1_ )->being_packed( resi ) || stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
				// If taking the relative complement of task1 in task2, then prevent repacking at all residues that are in task2 and not task1 only.
			} else if ( !invert_ && operator_ == "NOT" ) {
				if ( !(stored_tasks.get_task( task1_ )->being_packed( resi ) && !stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
				// If taking the relative complement of task2 in task1, then prevent repacking at all residues that are in task1 and not task2 only.
			} else if ( invert_ && operator_ == "NOT" ) {
				if ( (stored_tasks.get_task( task1_ )->being_packed( resi ) && !stored_tasks.get_task( task2_ )->being_packed( resi )) ) {
					task->nonconst_residue_task(resi).prevent_repacking();
				} else {
					select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
				}
			} else {
				utility_exit_with_message("Error: " + operator_ + " is not a valid option for the logic of the StoreCombinedStoredTaskMover.");
			}
		}
		TR << select_packable_pos << std::endl;
		// Store the new combined tasks in the cachebable data of the pose.
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
		}
		// If you haven't set overwrite to true and your task name already exists, fail. Otherwise, put the task you've made into the data cache.
		if ( overwrite_ || !stored_tasks.has_task(task_name_) ) {
			stored_tasks.set_task( task, task_name_ );
		} else {
			utility_exit_with_message("A stored task with the name " + task_name_ + " already exists; you must set overwrite flag to true to overwrite." );
		}
	}
}

void
StoreCombinedStoredTasksMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & /*data_map*/, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task1_ = tag->getOption< std::string >( "task1" ) ;
	task2_ = tag->getOption< std::string >( "task2" ) ;
	operator_ = tag->getOption< std::string >( "operator", "union" ) ;
	task_name_ = tag->getOption< std::string >( "task_name" );
	overwrite_ = tag->getOption< bool >( "overwrite", 0 );
	invert_ = tag->getOption< bool >( "invert", 0 );
}

// @brief Identification
std::string StoreCombinedStoredTasksMoverCreator::keyname() const { return StoreCombinedStoredTasksMoverCreator::mover_name(); }
std::string StoreCombinedStoredTasksMoverCreator::mover_name() { return "StoreCombinedStoredTasksMover"; }
std::string StoreCombinedStoredTasksMover::get_name() const { return "StoreCombinedStoredTasksMover"; }

protocols::moves::MoverOP
StoreCombinedStoredTasksMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StoreCombinedStoredTasksMover );
}

protocols::moves::MoverOP
StoreCombinedStoredTasksMover::clone() const {
	return protocols::moves::MoverOP( new StoreCombinedStoredTasksMover( *this ) );
}

protocols::moves::MoverOP
StoreCombinedStoredTasksMover::fresh_instance() const {
	return protocols::moves::MoverOP( new StoreCombinedStoredTasksMover );
}

} // task_operations
} // toolbox
} // protocols

