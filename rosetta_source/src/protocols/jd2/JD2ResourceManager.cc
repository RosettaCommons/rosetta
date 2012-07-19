// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManager.cc
/// @brief  The ResourceManager that is compatible with the JD2 JobDistributor
/// @author Andrew Leaver-Fay
/// @author Brian Weitzner
/// @author Matthew O'Meara

//unit headers
#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerCreator.hh>
#include <basic/resource_manager/types.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/FallbackConfigurationFactory.hh>
#include <basic/resource_manager/JobOptions.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Project Headers
#include <basic/options/option.hh>

// basic Headers
#include <basic/Tracer.hh>
#include <basic/resource_manager/types.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/ResourceLocatorFactory.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/ResourceOptionsFactory.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//C++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace jd2 {

using std::stringstream;
using std::string;
using basic::Tracer;
using basic::resource_manager::ResourceManager;
using basic::resource_manager::ResourceDescription;
using basic::resource_manager::ResourceOP;
using basic::resource_manager::JobOptionsOP;
using utility::excn::EXCN_Msg_Exception;

static Tracer TR("protocols.resource_manager.planner.JD2ResourceManager");

JD2ResourceManagerCreator::~JD2ResourceManagerCreator() {}

ResourceManager *
JD2ResourceManagerCreator::create_resource_manager() const
{
	return new JD2ResourceManager;
}


string
JD2ResourceManagerCreator::manager_type() const
{
	return "JD2ResourceManager";
}

void
JD2ResourceManager::clear()
{
	basic::resource_manager::LazyResourceManager::clear();
}


JD2ResourceManager::JD2ResourceManager() {
}

JD2ResourceManager::~JD2ResourceManager() {}

/// @details instantiate all the resource locators given in the input tags, and
/// put them into the base class.  Make sure no two resource locators share a common
/// name.
void JD2ResourceManager::read_resource_locators_tags( utility::tag::TagPtr tags )
{
	using basic::resource_manager::ResourceLocatorFactory;
	using basic::resource_manager::ResourceLocatorOP;
	using basic::resource_manager::LocatorTag;
	using utility::tag::Tag;
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	for ( Tag::tags_t::const_iterator
			tag_iter = tags->getTags().begin(),
			tag_iter_end = tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();

		/// 1. Make sure the resource locator has been given a tag.
		if ( ! (*tag_iter)->hasOption( "tag" ) ) {
			std::ostringstream err;
			err << "Unable to find a 'tag' for a ResourceLocator of type '" << tagname << "'\n";
			throw MsgException( err.str() );
		}
		LocatorTag locator_tag = (*tag_iter)->getOption< LocatorTag >( "tag" );

		// 2. Make sure no other resource_options object has been declared with this name
		if ( LazyResourceManager::has_resource_locator( locator_tag ) ) {
			std::ostringstream err;
			err << "Duplicated tag, '" << locator_tag <<"' assigned to a ResourceLocator object with type '";
			err << tagname << "'.\n";
			err << "Prevous ResourceLocator object with this tag was of type '";
			err << LazyResourceManager::find_resource_options( locator_tag )->type() << "'\n";
			throw MsgException( err.str() );
		}

		/// 3. Try to create this ResourceLocator object; the factory may throw.  Catch any thrown MsgException and
		/// append to its message the tagname and options_tag for the ResourceOptions being read.
		ResourceLocatorOP resource_options;
		try {
			resource_options = ResourceLocatorFactory::get_instance()->create_resource_locator( tagname, *tag_iter );
		} catch ( MsgException const & e ) {
			std::ostringstream err;
			err << e.msg() << "\n";
			err << "Exception thrown while trying to initialize a ResourceLocator of type '";
			err << tagname << "' with a tag of '" << locator_tag << "'\n";
			err << "from JD2ResourceManager::read_resource_locators_tags\n";
			throw MsgException( err.str() );
		}
		LazyResourceManager::add_resource_locator( locator_tag, resource_options );
	}

}

/// @details instantiate all the resource options and put them in the base class.
/// Make sure no two resource options are given the same name.
void JD2ResourceManager::read_resource_options_tags( utility::tag::TagPtr tags )
{
	using basic::resource_manager::ResourceOptionsFactory;
	using basic::resource_manager::ResourceOptionsOP;
	using basic::resource_manager::ResourceOptionsTag;
	using utility::tag::Tag;
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	for ( Tag::tags_t::const_iterator
			tag_iter = tags->getTags().begin(),
			tag_iter_end = tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();

		/// 1. Make sure this resource_options object has been declared.
		if ( ! (*tag_iter)->hasOption( "tag" ) ) {
			std::ostringstream err;
			err << "Unable to find a 'tag' for a ResourceOption of type '" << tagname << "'\n";
			throw MsgException( err.str() );
		}
		ResourceOptionsTag options_tag = (*tag_iter)->getOption< ResourceOptionsTag >( "tag" );

		// 2. Make sure no other resource_options object has been declared with this name
		if ( LazyResourceManager::has_resource_options( options_tag ) ) {
			std::ostringstream err;
			err << "Duplicated tag, '" << options_tag <<"' assigned to a ResourceOptions object with type '";
			err << tagname << "'.\n";
			err << "Prevous ResourceOptions object with this tag was of type '";
			err << LazyResourceManager::find_resource_options( options_tag )->type() << "'\n";
			throw MsgException( err.str() );
		}

		/// 3. Try to create this ResourceOptions object; the factory may throw.  Catch any thrown MsgException and
		/// append to its message the tagname and options_tag for the ResourceOptions being read.
		ResourceOptionsOP resource_options;
		try {
			resource_options = ResourceOptionsFactory::get_instance()->create_resource_options( tagname, *tag_iter );
		} catch ( MsgException const & e ) {
			std::ostringstream err;
			err << e.msg() << "\n";
			err << "Exception thrown while trying to initialize a ResourceOption of type '";
			err << tagname << "' with a tag of '" << options_tag << "'\n";
			err << "from JD2ResourceManager::read_resource_options_tags\n";
			throw MsgException( err.str() );
		}
		LazyResourceManager::add_resource_options( options_tag, resource_options );
	}

}

/// @details read through all the resources, and put them into the base class
/// for later instantiation.  Make sure each resource is named, that it is the
/// only resource that has been declared with this name, and that
void JD2ResourceManager::read_resources_tags( utility::tag::TagPtr tags )
{
	using basic::resource_manager::LoaderType;
	using basic::resource_manager::LocatorTag;
	using basic::resource_manager::LocatorID;
	using basic::resource_manager::ResourceConfiguration;
	using basic::resource_manager::ResourceOptionsTag;
	using basic::resource_manager::ResourceLoaderFactory;
	using basic::resource_manager::ResourceTag;
	using utility::tag::Tag;
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	for ( Tag::tags_t::const_iterator
			tag_iter = tags->getTags().begin(),
			tag_iter_end = tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		LoaderType loader_type = (*tag_iter)->getName();

		/// 1. Make sure this is an allowed resource type / loader type
		if ( ! ResourceLoaderFactory::get_instance()->has_resource_loader( loader_type ) ) {
			std::ostringstream err;
			err << "ResourceLoader type '" << loader_type << "' requested in JD2ResourceManager::read_resources_tags has not been registered with the ResourceLoaderFactory.  Available types include:.\n";
			std::list< std::string > loader_types = ResourceLoaderFactory::get_instance()->available_resource_loaders();
			for ( std::list< std::string >::const_iterator
					iter = loader_types.begin(), iter_end = loader_types.end(); iter != iter_end; ++iter ) {
				err << "  '" << *iter << std::endl;
			}
			throw MsgException( err.str() );
		}

		/// 2. Make sure this resource object has been given a tag.
		if ( ! (*tag_iter)->hasOption( "tag" ) ) {
			std::ostringstream err;
			err << "Unable to find a 'tag' for a Resource of type '" << loader_type << "'\n";
			throw MsgException( err.str() );
		}
		ResourceTag resource_tag = (*tag_iter)->getOption< ResourceTag >( "tag" );

		// 3. Make sure no other resource object has been declared with this name
		if ( LazyResourceManager::has_resource_configuration( resource_tag ) ) {
			std::ostringstream err;
			err << "Duplicated tag, '" << resource_tag << "' assigned to a Resource object with ResourceLoader type '";
			err << loader_type << "'.\n";
			throw MsgException( err.str() );
		}

		/// 4. Verify that, if it's been given a locator, that the ResourceLocator has previously been declared
		LocatorTag locator_tag;
		if ( (*tag_iter)->hasOption( "locator" ) ) {
			locator_tag = (*tag_iter)->getOption< LocatorTag >( "locator" );
			if ( ! LazyResourceManager::has_resource_locator( locator_tag ) ) {
				std::ostringstream err;
				err << "Resource '" << resource_tag << "' with LoaderType '" << loader_type << "' was given the tag ";
				err << "for a ResourceLocator, '" << locator_tag << "', which has not previously been declared.";
				throw MsgException( err.str() );
			}
		}

		/// 5. Verify that, if it's been given a resource options, that the ResourceOptions has been previously declared
		ResourceOptionsTag resource_options_tag;
		if ( (*tag_iter)->hasOption( "options" ) ) {
			resource_options_tag = (*tag_iter)->getOption< ResourceOptionsTag >( "options" );
			if ( ! LazyResourceManager::has_resource_options( resource_options_tag ) ) {
				std::ostringstream err;
				err << "Resource '" << resource_tag << "' with LoaderType '" << loader_type << "' was given the tag ";
				err << "for a ResourceLoaderOptions, '" << resource_options_tag << "', which has not previously been declared.";
				throw MsgException( err.str() );
			}
		}

		/// 6. Verify that it has been given a locatorID
		if ( ! (*tag_iter)->hasOption( "locatorID" ) ) {
			std::ostringstream err;
			err << "Resource '" << resource_tag << "' with LoaderType '" << loader_type << "' was not supplied with ";
			err << "an option 'locatorID', which is required" << std::endl;
			throw MsgException( err.str() );
		}
		LocatorID locator_id = (*tag_iter)->getOption< LocatorID >( "locatorID" );

		ResourceConfiguration resource_configuration;
		resource_configuration.resource_tag           = resource_tag;
		resource_configuration.locator_tag            = locator_tag;
		resource_configuration.locator_id             = locator_id;
		resource_configuration.loader_type            = resource_tag;
		resource_configuration.resource_options_tag   = resource_options_tag;

		LazyResourceManager::add_resource_configuration( resource_tag, resource_configuration );
	}

}


JD2ResourceManager *
JD2ResourceManager::get_jd2_resource_manager_instance(
){
	JD2ResourceManager * jd2_resource_manager(
		dynamic_cast< JD2ResourceManager * > ( ResourceManager::get_instance() ));
	if ( ! jd2_resource_manager ) {
		std::ostringstream err;
		err
			<< "Error trying to access the JD2ResourceManager "
			<< "from the JD2ResourceManagerJobInputer";
		EXCN_Msg_Exception( err.str() );
	}
	return jd2_resource_manager;
}

bool
JD2ResourceManager::has_resource_with_description(
	ResourceDescription const & resource_description
) {
	using basic::resource_manager::FallbackConfigurationFactory;
	using basic::resource_manager::FallbackConfigurationOP;
	using basic::resource_manager::ResourceTag;
	
	if ( has_resource_tag_by_job_tag(
		resource_description,
		JobDistributor::get_instance()->current_job()->input_tag()) )
	{
		return true;
	}
	
	if ( FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( resource_description ))
	{
		FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( resource_description );
		ResourceTag resource_tag = fallback->get_resource_tag_from_description( resource_description );
		add_resource_configuration( resource_tag, create_resource_configuration_from_fallback( fallback, resource_description ));
		add_resource_tag_by_job_tag( resource_description, JobDistributor::get_instance()->current_job()->input_tag(), resource_tag );
		return true;
	}
	return false;
}

ResourceOP
JD2ResourceManager::get_resource(
	ResourceDescription const & resource_description
) {
	return get_resource_by_job_tag(
		resource_description,
		JobDistributor::get_instance()->current_job()->input_tag());
}

bool
JD2ResourceManager::get_option(
	utility::options::BooleanOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< bool > const &
JD2ResourceManager::get_option(
	utility::options::BooleanVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::file::FileName const &
JD2ResourceManager::get_option(
	utility::options::FileOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< utility::file::FileName > const &
JD2ResourceManager::get_option(
	utility::options::FileVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

int
JD2ResourceManager::get_option(
	utility::options::IntegerOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< int > const &
JD2ResourceManager::get_option(
	utility::options::IntegerVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::file::PathName const &
JD2ResourceManager::get_option(
	utility::options::PathOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< utility::file::PathName > const &
JD2ResourceManager::get_option(
	utility::options::PathVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

platform::Real
JD2ResourceManager::get_option(
	utility::options::RealOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< platform::Real > const &
JD2ResourceManager::get_option(
	utility::options::RealVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

std::string const &
JD2ResourceManager::get_option(
	utility::options::StringOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

utility::vector1< std::string > const &
JD2ResourceManager::get_option(
	utility::options::StringVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( has_job_options( currjob.input_tag() ) ) {
		JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
		if(job_option->has_option(key)){
			return job_option->get_option(key);
		}
	}
	return basic::options::option[ key ]();
}

bool
JD2ResourceManager::has_option(
	utility::options::BooleanOptionKey key
) const {
	
	if ( ! has_job_options( JobDistributor::get_instance()->current_job()->input_tag() )) return false;
	JobOptionsOP job_option(
		get_job_options(
			JobDistributor::get_instance()->current_job()->input_tag()));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::BooleanVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::FileOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::FileVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::IntegerOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::IntegerVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::PathOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::PathVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::RealOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::RealVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::StringOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

bool
JD2ResourceManager::has_option(
	utility::options::StringVectorOptionKey key
) const {
	Job const &  currjob( * JobDistributor::get_instance()->current_job() );
	if ( ! has_job_options( currjob.input_tag() ) ) { return false; }
	JobOptionsOP job_option( get_job_options( currjob.input_tag() ));
	return job_option->has_option(key);
}

basic::resource_manager::ResourceConfiguration
JD2ResourceManager::create_resource_configuration_from_fallback(
	basic::resource_manager::FallbackConfigurationCOP fallback,
	ResourceDescription const & resource_description
)
{
	using basic::resource_manager::ResourceConfiguration;
	ResourceConfiguration resource_configuration;
	resource_configuration.resource_tag = fallback->get_resource_tag_from_description( resource_description );
	resource_configuration.locator_tag = fallback->get_locator_tag_from_description( resource_description );
	resource_configuration.locator_id = fallback->get_locator_id_from_description( resource_description );
	resource_configuration.loader_type = fallback->get_loader_type_from_description( resource_description );
	resource_configuration.resource_options_tag = fallback->get_resource_options_tag_from_description( resource_description );
	return resource_configuration;
}

} // namespace
} // namespace
