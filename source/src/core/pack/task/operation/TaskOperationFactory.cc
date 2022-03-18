// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/TaskOperationFactory.cc
/// @brief
/// @author ashworth

#include <core/pack/task/operation/TaskOperationFactory.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>
#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationManager.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/vector0.hh>

// Boost headers

namespace core {
namespace pack {
namespace task {
namespace operation {

static basic::Tracer TR( "core.pack.task.operation.TaskOperationFactory" );

TaskOperationFactory::~TaskOperationFactory()= default;

/// @brief the default TaskOperations are now initialized in core/init/init.cc via the registrator/creator scheme
TaskOperationFactory::TaskOperationFactory() = default;

void
TaskOperationFactory::factory_register( TaskOperationCreatorOP creator )
{
	if ( task_operation_creator_map_.find( creator->keyname() ) != task_operation_creator_map_.end() ) {
		utility_exit_with_message( "Factory Name Conflict: Two or more TaskOperationCreators registered with the name " + creator->keyname() );
	}
	add_creator( creator );
}

/// @brief add a TaskOperation prototype creator
void
TaskOperationFactory::add_creator( TaskOperationCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	task_operation_creator_map_[ creator->keyname() ] = creator;
}

bool TaskOperationFactory::has_type( std::string const & type ) const
{
	return ( task_operation_creator_map_.find( type ) != task_operation_creator_map_.end() );
}

/// @brief adds a ResLvlTaskOperation prototype creator to the child ResLvlTaskOperationFactory
void
TaskOperationFactory::add_creator( ResLvlTaskOperationCreatorOP creator )
{
	ResLvlTaskOperationFactory::get_instance()->add_creator( creator );
}

/// @brief adds a ResFilter prototype creator to the child ResFilterFactory
void
TaskOperationFactory::add_creator( ResFilterCreatorOP creator )
{
	ResFilterFactory::get_instance()->add_creator( creator );
}

/// @brief return new TaskOperation by key lookup in task_operation_creator_map_ (new TaskOperation parses Tag if provided)
/*!
Example Tag syntax for parser as of Summer 2009

<ReadResfile name=rrf filename=myresfile/>

or

<OperateOnCertainResidues name=PROTEINnopack>
<PreventRepackingRLT/>
<ResidueHasProperty property=PROTEIN/>
</OperateOnCertainResidues>

*/

/// @brief Get the XML schema for a given TaskOperation.
/// @details Throws an error if the TaskOperation is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
TaskOperationFactory::provide_xml_schema(
	std::string const &task_operation_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	auto iter( task_operation_creator_map_.find( task_operation_name ) );
	if ( iter != task_operation_creator_map_.end() ) {
		iter->second->provide_xml_schema( xsd );
	} else {
		TR << "Available options: ";
		for ( const auto & to_iter : task_operation_creator_map_ ) {
			TR << to_iter.first << ", ";
		}
		TR << std::endl;
		utility_exit_with_message( task_operation_name + " is not known to the TaskOperationFactory. Was it registered in the appropriate initialization files (src/protocols/init/init.TaskOperationCreators.ihh and src/protocols/init/init.TaskOperationRegistrators.ihh)?" );
	}
}

/// @brief Get a human-readable listing of the citations for a given filter, by filter name.
/// @details Returns an empty string if there are no citations.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
TaskOperationFactory::get_citation_humanreadable(
	std::string const & taskop_name
) const {
	using namespace basic::citation_manager;
	CitationCollectionList citations;

	auto const iter( task_operation_creator_map_.find( taskop_name ) );
	runtime_assert_string_msg( iter != task_operation_creator_map_.end(), "Error in TaskOperationFactory::get_citation_humanreadable(): Could not parse \"" + taskop_name + "\"!" );

	TaskOperationOP task_operation( iter->second->create_task_operation() );
	runtime_assert_string_msg( task_operation != nullptr, "Error in TaskOperationFactory::get_citation_humanreadable(): Could not instantiate " + taskop_name + "!" );
	task_operation->provide_citation_info(citations);
	if ( citations.empty() ) return "";
	std::ostringstream ss;
	ss << "References and author information for the " << taskop_name << " task operation:" << std::endl;
	ss << std::endl;
	basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info_from_list_to_stream( citations, ss );
	return ss.str();
}

TaskOperationOP
TaskOperationFactory::newTaskOperation(
	std::string const & type,
	basic::datacache::DataMap & datamap,
	TagCOP tag /* = boost::shared_ptr< Tag >() */
) const
{
	auto iter( task_operation_creator_map_.find( type ) );
	if ( iter != task_operation_creator_map_.end() ) {
		TaskOperationOP task_operation( iter->second->create_task_operation() );
		// parse tag if tag pointer is pointing to one
		if ( tag.get() != nullptr ) task_operation->parse_tag( tag, datamap );

		// Register this task operation with the citation manager:
		basic::citation_manager::CitationCollectionList citations;
		task_operation->provide_citation_info( citations );
		basic::citation_manager::CitationManager::get_instance()->add_citations( citations );

		return task_operation;
	} else {
		TR<<"Available options: ";
		for ( const auto & to_iter : task_operation_creator_map_ ) {
			TR<<to_iter.first<<", ";
		}
		TR<<std::endl;
		utility_exit_with_message( type + " is not known to the TaskOperationFactory. Was it registered in the appropriate initialization files (src/protocols/init/init.TaskOperationCreators.ihh and src/protocols/init/init.TaskOperationRegistrators.ihh)?" );
		return nullptr;
	}
}

/// @brief recurse tag file to find TASKOPERATIONS definitions
void
TaskOperationFactory::newTaskOperations( TaskOperationOPs & tops, basic::datacache::DataMap & datamap, TagCOP tag ) const
{
	using TagCOPs = utility::vector0<TagCOP>;
	TR.Trace << "Tag name " << tag->getName();
	if ( tag->getTags().empty() ) { TR.Trace << " (empty)" << std::endl; return; }
	else TR.Trace << std::endl;
	TagCOPs const subtags( tag->getTags() );
	if ( tag->getName() == "TASKOPERATIONS" ) {
		for ( const auto & subtag : subtags ) {
			std::string const type( subtag->getName() );
			TaskOperationOP new_to = newTaskOperation( type, datamap, subtag );
			runtime_assert( new_to != nullptr );
			tops.push_back( new_to );
			TR << "Created and parsed anonymous TaskOperation of type " << type << std::endl;
		}
	}
	// recurse
	for ( const auto & subtag : subtags ) {
		newTaskOperations( tops, datamap, subtag );
	}
}

void
TaskOperationFactory::newTaskOperations( TaskOperationOPs & tops, basic::datacache::DataMap & datamap, std::string const & tagfilename ) const
{
	utility::io::izstream fin;
	fin.open( tagfilename.c_str() );
	runtime_assert( fin.good() );
	TagCOP tag = utility::tag::Tag::create(fin);
	fin.close();
	TR << "TaskOperationFactory parsing " << tagfilename << " to create TaskOperations:" << std::endl;
	TR << tag << std::endl;
	newTaskOperations( tops, datamap, tag );
}

void TaskOperationFactory::define_task_op_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			task_operation_creator_map_,
			task_operation_xml_schema_group_name(),
			& complex_type_name_for_task_op,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for TaskOperations from TaskOperationFactory\n" + e.msg() );
	}

}

TaskOperationFactory::TaskOperationCreatorMap const & TaskOperationFactory::creator_map() const
{
	return task_operation_creator_map_;
}

std::string TaskOperationFactory::task_operation_xml_schema_group_name() {
	return "task_operation";
}

TaskOperationFactory::TagCOP TaskOperationFactory::MakeTagCOP()
{
	return utility::pointer::make_shared< utility::tag::Tag >();
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
