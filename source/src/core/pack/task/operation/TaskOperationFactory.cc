// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperationFactory.cc
/// @brief
/// @author ashworth

#include <core/pack/task/operation/TaskOperationFactory.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/ResFilterCreator.hh>
#include <core/pack/task/operation/ResFilterFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <utility/vector0.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace pack {
namespace task {
namespace operation {

static THREAD_LOCAL basic::Tracer TR( "core.pack.task.operation.TaskOperationFactory" );

// special singleton functions
// initialize
#if defined MULTI_THREADED && defined CXX11
std::atomic< TaskOperationFactory * > TaskOperationFactory::instance_( 0 );
#else
TaskOperationFactory * TaskOperationFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex TaskOperationFactory::singleton_mutex_;

std::mutex & TaskOperationFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
TaskOperationFactory * TaskOperationFactory::get_instance()
{
	boost::function< TaskOperationFactory * () > creator = boost::bind( &TaskOperationFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

TaskOperationFactory *
TaskOperationFactory::create_singleton_instance()
{
	return new TaskOperationFactory;
}

TaskOperationFactory::~TaskOperationFactory(){}

/// @brief the default TaskOperations are now initialized in core/init/init.cc via the registrator/creator scheme
TaskOperationFactory::TaskOperationFactory() {}

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
	runtime_assert( creator != 0 );
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
TaskOperationOP
TaskOperationFactory::newTaskOperation(
	std::string const & type,
	basic::datacache::DataMap & datamap,
	TagCOP tag /* = boost::shared_ptr< Tag >() */
) const
{
	TaskOperationCreatorMap::const_iterator iter( task_operation_creator_map_.find( type ) );
	if ( iter != task_operation_creator_map_.end() ) {
		TaskOperationOP task_operation( iter->second->create_task_operation() );
		// parse tag if tag pointer is pointing to one
		if ( tag.get() != NULL ) task_operation->parse_tag( tag, datamap );
		return task_operation;
	} else {
		TR<<"Available options: ";
		for ( TaskOperationCreatorMap::const_iterator to_iter = task_operation_creator_map_.begin(); to_iter != task_operation_creator_map_.end(); ++to_iter ) {
			TR<<to_iter->first<<", ";
		}
		TR<<std::endl;
		utility_exit_with_message( type + " is not known to the TaskOperationFactory. Was its taskOperationCreator class registered at initialization?" );
		return NULL;
	}
}

/// @brief recurse tag file to find TASKOPERATIONS definitions
void
TaskOperationFactory::newTaskOperations( TaskOperationOPs & tops, basic::datacache::DataMap & datamap, TagCOP tag ) const
{
	typedef utility::vector0< TagCOP > TagCOPs;
	TR.Trace << "Tag name " << tag->getName();
	if ( tag->getTags().empty() ) { TR.Trace << " (empty)" << std::endl; return; }
	else TR.Trace << std::endl;
	TagCOPs const subtags( tag->getTags() );
	if ( tag->getName() == "TASKOPERATIONS" ) {
		for ( TagCOPs::const_iterator tp( subtags.begin() ), tp_e( subtags.end() ); tp != tp_e; ++tp ) {
			std::string const type( (*tp)->getName() );
			TaskOperationOP new_to = newTaskOperation( type, datamap, *tp );
			runtime_assert( new_to != 0 );
			tops.push_back( new_to );
			TR << "Created and parsed anonymous TaskOperation of type " << type << std::endl;
		}
	}
	// recurse
	for ( TagCOPs::const_iterator tp( subtags.begin() ), tp_e( subtags.end() ); tp != tp_e; ++tp ) {
		newTaskOperations( tops, datamap, *tp );
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

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
