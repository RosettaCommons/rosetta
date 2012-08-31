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

// Numeric headers
#include <numeric/random/random.hh>

// basic Headers
#include <basic/Tracer.hh>
#include <basic/resource_manager/types.hh>
#include <basic/resource_manager/ResourceLoader.hh>
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
#include <utility/string_util.hh>

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
using basic::resource_manager::LoaderType;
using basic::resource_manager::LocatorTag;
using basic::resource_manager::LocatorID;
using basic::resource_manager::ResourceConfiguration;
using basic::resource_manager::ResourceOptionsTag;
using basic::resource_manager::ResourceLoaderFactory;
using basic::resource_manager::ResourceTag;
using utility::tag::Tag;

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

		// 2. Make sure no other resource locator object has been declared with this name
		if ( LazyResourceManager::has_resource_locator( locator_tag ) ) {
			std::ostringstream err;
			err << "Duplicated tag, '" << locator_tag <<"' assigned to a ResourceLocator object with type '";
			err << tagname << "'.\n";
			err << "Prevous ResourceLocator object with this tag was of type '";
			err << LazyResourceManager::find_resource_locator( locator_tag )->type() << "'\n";
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

///@detail Create the loader type from the name of the resource tag
LoaderType
JD2ResourceManager::read_resource_loader_type_item(
	utility::tag::TagPtr tag
) {
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	LoaderType loader_type = tag->getName();

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
	return loader_type;
}

///@details make sure the resource object has been given a tag and
///that no other resource object has been delecared with the same name.
ResourceTag
JD2ResourceManager::read_resource_tag_item(
	utility::tag::TagPtr tag,
	LoaderType const & loader_type,
	LocatorID const & locator_id
) {
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	/// 2. Set the resource tag to be the locator id unless it has been
	/// explcitely provided.
	ResourceTag resource_tag;
	if ( ! tag->hasOption( "tag" ) ) {
		resource_tag = locator_id;
	} else {
		resource_tag = tag->getOption< ResourceTag >( "tag" );
	}

	// 3. Make sure no other resource object has been declared with this name
	if ( LazyResourceManager::has_resource_configuration( resource_tag ) ) {
		std::ostringstream err;
		err << "Duplicated tag, '" << resource_tag << "' assigned to a Resource object with ResourceLoader type '";
		err << loader_type << "'.\n";
		throw MsgException( err.str() );
	}
	return resource_tag;
}

///@brief find the LocatorTag item from the resource tag. Based on the
///LocatorTag fill the LocatorID.
///   locator item (&string):
///     FileSystemResourceLocator (default)
///       - file item is interchangable with the locatorID item
///   locatorID item (&string):
LocatorTag
JD2ResourceManager::read_resource_locator_items(
	utility::tag::TagPtr tag,
	LoaderType const & loader_type,
	LocatorID & locator_id
) {
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	/// 4. Verify that, if it's been given a locator, that the ResourceLocator has previously been declared
	LocatorTag locator_tag;
	if ( tag->hasOption( "locator" ) ) {
		locator_tag = tag->getOption< LocatorTag >( "locator" );
		if ( ! LazyResourceManager::has_resource_locator( locator_tag ) ) {
			std::ostringstream err;
			err
				<< "Resource subtag"
				<< " with LoaderType '" << loader_type << "'"
				<< " was given the ResourceLocator tag '" << locator_tag << "',"
				<< " which has not previously been declared.";
			throw MsgException( err.str() );
		}
		if( tag->hasOption( "file" ) && locator_tag != "FileSystemResourceLocator" ){
			std::ostringstream err;
			err
				<< "Resource subtag "
				<< "with LoaderType '" << loader_type << "'"
				<< " has locator='" << locator_tag << "' and "
				<< " file='" << tag->getOption< LocatorID > ("file") << "'."
				<< " But specifying a file is only compatible with the"
				<< " FileSystemResourceLocator." << std::endl;
			throw MsgException(err.str());
		}
	} else {
		locator_tag = "";
	}

	if(tag->hasOption( "file" )){
		if(tag->hasOption("locatorID")) {
			std::ostringstream err;
			err
				<< "Resource subtag"
				<< " with LoaderType '" << loader_type << "'"
				<< " has both"
				<< " file='" << tag->getOption< LocatorID >("file") << "' and"
				<< " locatorID='" << tag->getOption< LocatorID >("locatorID") << "'"
				<< " but it is only allowed to have one." << std::endl;
			throw MsgException(err.str());
		}
		locator_id = tag->getOption< LocatorID >("file");
	} else if( tag->hasOption("locatorID")){
		locator_id = tag->getOption< LocatorID >( "locatorID" );
	} else {
		std::ostringstream err;
		err
			<< "Resource subtag"
			<< " with LoaderType '" << loader_type << "' was not supplied"
			<< " with a locatorID tag, which is required." << std::endl;
		throw MsgException( err.str() );
	}

	return locator_tag;
}


ResourceOptionsTag
JD2ResourceManager::read_resource_options_tag_item(
	utility::tag::TagPtr tag,
	LoaderType const & loader_type,
	ResourceTag const & resource_tag
) {
	typedef utility::excn::EXCN_Msg_Exception MsgException;

	/// 5. Verify that, if it's been given a resource options, that the
	/// ResourceOptions has been previously declared
	ResourceOptionsTag resource_options_tag;
	if ( tag->hasOption( "options" ) ) {
		resource_options_tag = tag->getOption< ResourceOptionsTag >( "options" );
		if ( ! LazyResourceManager::has_resource_options( resource_options_tag ) ) {
			std::ostringstream err;
			err << "Resource '" << resource_tag << "' with LoaderType '" << loader_type << "' was given the tag ";
			err << "for a ResourceLoaderOptions, '" << resource_options_tag << "', which has not previously been declared.";
			throw MsgException( err.str() );
		}
	}
	return resource_options_tag;
}

/// @details read through all the resources, and put them into the base class
/// for later instantiation.  Make sure each resource is named, that it is the
/// only resource that has been declared with this name, and that
void JD2ResourceManager::read_resources_tags( utility::tag::TagPtr tags )
{

	for ( Tag::tags_t::const_iterator
			tag_iter = tags->getTags().begin(),
			tag_iter_end = tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {

		LoaderType loader_type(read_resource_loader_type_item(*tag_iter));

		LocatorID locator_id;
		LocatorTag locator_tag(
			read_resource_locator_items(
				*tag_iter, loader_type, locator_id));

		ResourceTag resource_tag(
			read_resource_tag_item(*tag_iter, loader_type, locator_id));


		ResourceOptionsTag resource_options_tag(
			read_resource_options_tag_item(*tag_iter, loader_type, resource_tag));

		ResourceConfiguration resource_configuration;
		resource_configuration.loader_type            = loader_type;
		resource_configuration.resource_tag           = resource_tag;
		resource_configuration.locator_tag            = locator_tag;
		resource_configuration.locator_id             = locator_id;
		resource_configuration.resource_options_tag   = resource_options_tag;

		LazyResourceManager::add_resource_configuration(
			resource_tag, resource_configuration );
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
		JobDistributor::get_instance()->current_job()->input_tag()) ) {
		return true;
	}

	/// check, have we already created a fallback resource description?
	if ( fallback_resource_descriptions_created_.find( resource_description )
			!= fallback_resource_descriptions_created_.end() ) {
		return true;
	}

	if ( FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( resource_description )) {
		FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( resource_description );
		return fallback->fallback_specified( resource_description );
	}
	return false;
}

ResourceOP
JD2ResourceManager::get_resource(
	ResourceDescription const & resource_description
) {
	std::string const & jobtag = JobDistributor::get_instance()->current_job()->input_tag();

	if ( has_resource_tag_by_job_tag( resource_description, jobtag ) ) {
		return get_resource_by_job_tag( resource_description, jobtag);
	}

	std::map< std::string, std::string >::const_iterator fallbackname_iterator
		= fallback_resource_descriptions_created_.find( resource_description );

	if ( fallbackname_iterator !=  fallback_resource_descriptions_created_.end() ) {
		return find_resource( fallbackname_iterator->second );
	} else {
		using basic::resource_manager::FallbackConfigurationFactory;
		using basic::resource_manager::FallbackConfigurationOP;

        if ( FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( resource_description )) {
            FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( resource_description );
            if ( fallback->fallback_specified( resource_description ) ) {
                ResourceOP fallbackresource = create_resource_from_fallback( fallback, resource_description );
                // now make sure that the resource is saved for later; create a pheax uuid for this
                std::string fallbackname = "fallback_" + resource_description;
                for ( core::Size ii = 1; ii <= 10000; ++ii ) {
                    if ( ! has_resource_configuration( fallbackname )) break;
                    fallbackname = "fallback_" + resource_description + "_" + utility::to_string( numeric::random::random_range(1,4000000) );
                }
                if ( has_resource_configuration( fallbackname ) ) {
                    throw utility::excn::EXCN_Msg_Exception( "Could not name the fallback resource after 10000 attempts.  Try not to name your resources"
                        " beginning with the prefix 'fallback_'." );
                }
                basic::resource_manager::ResourceConfiguration fake_configuration;
                fake_configuration.resource_tag = fallbackname;
                add_resource_configuration( fallbackname, fake_configuration );
                fallback_resource_descriptions_created_[ resource_description ] = fallbackname;
                add_resource( fallbackname, fallbackresource );
                return fallbackresource;
            }
        }

		std::ostringstream errmsg;
		errmsg << "JD2ResourceManager does not have a resource "
			"corresponding to the resource description '" + resource_description + ".' for job '" +
			jobtag + "'.\n";
		errmsg << "Resources may be specified on the command line or through the JD2ResourceManagerJobInputter resource-declaration file.\n";

		if ( FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( resource_description ) ) {
			// i.e. this is a valid resource description, but the command line options required for this resource description
			// were not provided.
			errmsg << "The FallbackConfiguration for this resource description gives this error:\n";
			FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( resource_description );
			errmsg << fallback->could_not_create_resource_error_message( resource_description ) << "\n";
		} else {
			errmsg << "This resource description does not have a FallbackConfiguration defined.\n";
		}

		errmsg << "Thrown from JD2ResourceManager::get_resource\n";
		throw  utility::excn::EXCN_Msg_Exception( errmsg.str() );
	}
	return 0; // appease compiler
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

basic::resource_manager::ResourceOP
JD2ResourceManager::create_resource_from_fallback(
	basic::resource_manager::FallbackConfigurationCOP fallback,
	ResourceDescription const & resource_description
)
{
	using namespace basic::resource_manager;

	/// APL: Fix this.  For now, always use the FileSystemResourceLocator
	ResourceLocatorOP locator(find_resource_locator(""));

	ResourceStreamOP stream( locator->locate_resource_stream(
		fallback->get_locator_id(resource_description) ));

	ResourceLoaderOP loader(
		ResourceLoaderFactory::get_instance()->create_resource_loader(
			fallback->get_resource_loader( resource_description )));

	ResourceOptionsOP resource_options = fallback->get_resource_options( resource_description );
	if ( resource_options == 0 ) {
		resource_options = loader->default_options();
	}

	return loader->create_resource( *resource_options, fallback->get_locator_id( resource_description ), stream->stream());
}

} // namespace
} // namespace
