// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/TaskOperationLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace parser {

static thread_local basic::Tracer TR( "protocols.jd2.parser.TaskOperationLoader" );

TaskOperationLoader::TaskOperationLoader() {}
TaskOperationLoader::~TaskOperationLoader() {}

void TaskOperationLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace core::pack::task::operation;

	BOOST_FOREACH ( utility::tag::TagCOP tag, tag->getTags() ) {
		std::string const type( tag->getName() );
		if ( ! tag->hasOption("name") ) {
			utility_exit_with_message( "Can't create unnamed TaskOperation (type: " + type + ")" );
		}
		std::string const name( tag->getOption<std::string>("name") );
		if ( data.has( "task_operations", name ) ) {
			TR.Error << "Error TaskOperation of name \"" << name
				<< "\" (with type " << type << ") already exists. \n" << tag << std::endl;
			utility_exit_with_message("Duplicate definition of TaskOperation with name " + name);
		}
		TaskOperationOP new_t_o( TaskOperationFactory::get_instance()->newTaskOperation( type, data, tag ) );
		runtime_assert( new_t_o != 0 );
		data.add("task_operations", name, new_t_o );
		TR << "Defined TaskOperation named \"" << name << "\" of type " << type << std::endl;
	}
	TR.flush();
}

DataLoaderOP
TaskOperationLoaderCreator::create_loader() const { return DataLoaderOP( new TaskOperationLoader ); }

std::string
TaskOperationLoaderCreator::keyname() const { return "TASKOPERATIONS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
