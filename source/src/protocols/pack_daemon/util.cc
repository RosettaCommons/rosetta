// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/util.cc
/// @brief  Utility functions for working with the multistate-design classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/pack_daemon/util.hh>

// Package headers
#include <protocols/pack_daemon/EntityCorrespondence.hh>

// Core headers
#include <core/io/pdb/pose_io.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>

// Protocols headers
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/multistate_design/util.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace pack_daemon {

/// @throws Throws a utility::excn::EXCN_Msg_Exception if the input entity resfile does not
/// does not conform to the entity resfile format.
void
create_entity_resfile_contents(
	std::istream & resfile,
	std::string const & resfile_name,
	core::pack::task::ResfileContentsOP & entity_resfile_contents,
	core::pack::task::PackerTaskOP & entity_task,
	core::Size & num_entities
)
{
	using namespace core::pose;
	using namespace core::pack::task;

	if ( (resfile >> num_entities ).fail() ) {
		throw utility::excn::EXCN_Msg_Exception( "Error reading the number of entities from entity resfile " + resfile_name );
	}

	if ( num_entities == 0 ) {
		throw utility::excn::EXCN_Msg_Exception( "The number of entities in the input entity resfile must be greater than zero; error in " + resfile_name );
	}

	Pose ala_pose;
	std::string ala_seq( num_entities, 'A' );
	make_pose_from_sequence( ala_pose, ala_seq, core::chemical::FA_STANDARD );
	entity_task = TaskFactory::create_packer_task( ala_pose );

	entity_resfile_contents = new ResfileContents( ala_pose, resfile );

	/// apply the resfile operations to the entity_task_ for later error checking
	for ( core::Size ii = 1; ii <= entity_task->total_residue(); ++ii ) {

		std::list< ResfileCommandCOP > const & ii_command_list(
			entity_resfile_contents->specialized_commands_exist_for_residue( ii ) ?
			entity_resfile_contents->commands_for_residue( ii )
			: entity_resfile_contents->default_commands() );

		for ( std::list< ResfileCommandCOP >::const_iterator
				iter = ii_command_list.begin(), iter_end = ii_command_list.end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_action( *entity_task, ii );
		}
	}
}


void
initialize_task_from_entity_resfile_and_secondary_resfile(
	core::pose::Pose const & pose,
	EntityCorrespondenceCOP ec,
	core::pack::task::ResfileContents const & entity_resfile_contents,
	core::pack::task::ResfileContents const & secondary_resfile_contents,
	core::pack::task::PackerTaskOP task
)
{
	using namespace core::pose;
	using namespace core::pack::task;


	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		core::Size ii_entity = ec->entity_for_residue( ii );
		std::list< ResfileCommandCOP > const & ii_command_list(
			ii_entity != 0 ? (
				entity_resfile_contents.specialized_commands_exist_for_residue( ii_entity ) ?
				entity_resfile_contents.commands_for_residue( ii_entity ) :
				entity_resfile_contents.default_commands()
			) :
			(
				secondary_resfile_contents.specialized_commands_exist_for_residue( ii ) ?
				secondary_resfile_contents.commands_for_residue( ii ) :
				secondary_resfile_contents.default_commands() ));

		for ( std::list< ResfileCommandCOP >::const_iterator
				iter = ii_command_list.begin(), iter_end = ii_command_list.end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_action( *task, ii );
		}
	}
}

/// @details Create the list of available entity elements from the set of
/// available residue types in the "entity_task", the PackerTask that has
/// been created from an entity resfile, which defines a hypothetical pose
/// and places restrictions on the amino acids this pose can adopt.
/// These restrictions in turn must be translated into the available
/// search space that the genetic algorithm can explore.
void
initialize_ga_randomizer_from_entity_task(
	protocols::genetic_algorithm::PositionSpecificRandomizerOP rand,
	core::pack::task::PackerTaskOP entity_task
)
{
	for ( core::Size ii = 1; ii <= entity_task->total_residue(); ++ii ) {
		genetic_algorithm::EntityElements ii_elements =
			protocols::multistate_design::list_amino_acid_options(
			ii, entity_task->residue_task( ii ) );
		rand->append_choices( ii_elements );
	}
}

}
}

