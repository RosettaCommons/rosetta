// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/resource_manager/LazyResourceManager.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/ResourceLocatorFactory.hh>
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

//C++ headers
#include <string>
#include <sstream>
#include <iomanip>
#include <list>

namespace basic {
namespace resource_manager {

using std::stringstream;
using std::string;
using std::endl;
using std::setw;
using basic::Tracer;
using utility::vector1;

static Tracer TR("basic.resource_manager.LazyResourceManager");

void
ResourceConfiguration::show(
	std::ostream & out
) const {
	out
		<< "ResourceConfiguationKey Value" << endl
		<< "ResourceTag             " << resource_tag << endl
		<< "LocatorTag              " << locator_tag << endl
		<< "LocatorID               " << locator_id << endl
		<< "LoaderType              " << loader_type << endl
		<< "ResourceOptionsTag      " << resource_options_tag << endl;
}

std::ostream &
operator<<(
	std::ostream & out,
	ResourceConfiguration const & resource_configuration
) {
	resource_configuration.show(out);
	return out;
}


LazyResourceManager::LazyResourceManager() {
	add_default_resource_locator();
}

LazyResourceManager::~LazyResourceManager() = default;

void
LazyResourceManager::clear()
{
	resource_tags_.clear();
	resource_tag_lists_.clear();
	job_options_.clear();
	resource_configurations_.clear();
	resource_locators_.clear();
	resource_options_.clear();

	add_default_resource_locator();

}

void
LazyResourceManager::add_default_resource_locator()
{
	resource_locators_[""] = ResourceLocatorFactory::get_instance()->create_resource_locator( "FileSystemResourceLocator", "default_locator", nullptr );
	resource_locators_["NULL"] = ResourceLocatorFactory::get_instance()->create_resource_locator( "NullResourceLocator", "default_locator", nullptr );
}


///////////////////////////////////////////////////
///// ResourceDescription + JobTag interface //////
///////////////////////////////////////////////////

ResourceOP
LazyResourceManager::create_resource_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag
) const {
	ResourceTag resource_tag(find_resource_tag_by_job_tag(
		resource_description, job_tag));
	return create_resource(resource_tag);
}

void
LazyResourceManager::add_resource_tag_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag,
	ResourceTag const & resource_tag
)
{
	resource_tags_[make_pair(resource_description, job_tag)] = resource_tag;

	auto job_set_it(incomplete_job_sets_.find(resource_tag));
	if ( job_set_it == incomplete_job_sets_.end() ) {
		std::set<JobTag> new_set;
		new_set.insert(job_tag);
		incomplete_job_sets_[resource_tag]= new_set;
	} else {
		job_set_it->second.insert(job_tag);
	}

	auto resource_list_it(resource_tag_lists_.find(job_tag));
	if ( resource_list_it == resource_tag_lists_.end() ) {
		std::list<ResourceTag> new_list;
		new_list.push_back(resource_tag);
		resource_tag_lists_[job_tag] = new_list;
	} else {
		resource_list_it->second.push_back(resource_tag);
	}

}

bool
LazyResourceManager::has_resource_tag_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag
) const {
	auto resource_tag(
		resource_tags_.find(make_pair(resource_description, job_tag)));
	return resource_tag != resource_tags_.end();
}

ResourceTag
LazyResourceManager::find_resource_tag_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag
) const {
	auto resource_tag(
		resource_tags_.find(make_pair(resource_description, job_tag)));
	if ( resource_tag == resource_tags_.end() ) {
		stringstream err_msg;
		err_msg
			<< "Unable to find resource tag for the resource description '"
			<< resource_description << "' "
			<< "and job tag '" << job_tag << "'" << endl;
		utility_exit_with_message(err_msg.str());
	}
	return resource_tag->second;
}

ResourceOP
LazyResourceManager::get_resource_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag
) {
	ResourceTag resource_tag(find_resource_tag_by_job_tag(
		resource_description, job_tag));
	return find_resource(resource_tag);
}

std::list<ResourceTag>
LazyResourceManager::get_resource_tags_for_job_tag(
	JobTag const & job_tag
) const{
	auto resource_list_it(resource_tag_lists_.find(job_tag));
	if ( resource_list_it == resource_tag_lists_.end() ) {
		stringstream err_msg;
		err_msg
			<< "Unable to find job tag '" << job_tag << "'" << std::endl;
		utility_exit_with_message(err_msg.str());
	}
	return resource_list_it->second;
}


platform::Size
LazyResourceManager::get_count_of_jobs_associated_with_resource_tag(
	ResourceTag const & resource_tag) const
{
	auto job_set_it(incomplete_job_sets_.find(resource_tag));
	if ( job_set_it == incomplete_job_sets_.end() ) {
		stringstream err_msg;
		err_msg
			<< "Unable to find resource tag " << resource_tag <<std::endl;
		utility_exit_with_message(err_msg.str());
	}
	return job_set_it->second.size();
}

void
LazyResourceManager::mark_job_tag_as_complete(
	JobTag const & job_tag)
{


	std::list<ResourceTag> tag_list( get_resource_tags_for_job_tag(job_tag) );
	auto tag_list_it(tag_list.begin());
	for ( ; tag_list_it != tag_list.end(); ++tag_list_it ) {
		auto job_set_it(incomplete_job_sets_.find(*tag_list_it));
		auto job_tag_it(job_set_it->second.find(job_tag));
		if ( job_tag_it != job_set_it->second.end() ) {
			job_set_it->second.erase(job_tag_it);
		}
	}
}

void
LazyResourceManager::free_resource_by_tag(
	ResourceTag const & resource_tag
) {
	free_resource(resource_tag);
}

void
LazyResourceManager::free_resource_by_job_tag(
	ResourceDescription const & resource_description,
	JobTag const & job_tag
) {
	ResourceTag resource_tag(find_resource_tag_by_job_tag(
		resource_description, job_tag));
	free_resource(resource_tag);
}

///////////////////////////////////////////////////////////
///// Options ResourceDescription + JobTag interface //////
///////////////////////////////////////////////////////////

void
LazyResourceManager::add_job_options(
	JobTag const & job_tag,
	JobOptionsOP job_options
) {
	job_options_[job_tag] = job_options;
}

bool
LazyResourceManager::has_job_options(
	JobTag const & job_tag
) const {
	auto job_options(
		job_options_.find(job_tag));
	return job_options != job_options_.end();
}

JobOptionsOP
LazyResourceManager::get_job_options(
	JobTag const & job_tag
) const {
	auto job_options(
		job_options_.find(job_tag));
	if ( job_options == job_options_.end() ) {
		stringstream err_msg;
		err_msg
			<< "Unable to find job options for the resource description '"
			<< "' for the job with tag '" << job_tag << "'" << endl;
		utility_exit_with_message(err_msg.str());
	}
	return job_options->second;
}

//////////////////////////////////////////////
///// Resource Configuration interface  //////
//////////////////////////////////////////////

void
LazyResourceManager::add_resource_configuration(
	ResourceTag const & resource_tag,
	ResourceConfiguration const & resource_configuration
) {
	ResourceConfigurationMap::const_iterator config( resource_configurations_.find( resource_tag ));
	if ( config != resource_configurations_.end() ) {
		throw utility::excn::EXCN_Msg_Exception("Attempting to add multiple resource configurations with the resource tag '" + resource_tag + "'.");
	}
	resource_configurations_[ resource_tag ] = resource_configuration;
}

void
LazyResourceManager::add_resource_locator(
	LocatorTag const & locator_tag,
	ResourceLocatorOP resource_locator
) {
	ResourceLocatorsMap::const_iterator locator( resource_locators_.find( locator_tag ));
	if ( locator != resource_locators_.end() ) {
		throw utility::excn::EXCN_Msg_Exception("Attempting to add multiple resource locators with the locator tag '" + locator_tag + "'.");
	}
	resource_locators_[ locator_tag ] = resource_locator;
}


void
LazyResourceManager::add_resource_options(
	ResourceOptionsTag const & resource_options_tag,
	ResourceOptionsOP resource_options
) {
	ResourceOptionsMap::const_iterator options( resource_options_.find( resource_options_tag ));
	if ( options != resource_options_.end() ) {
		throw utility::excn::EXCN_Msg_Exception("Attempting to add multiple resource options with the resource optiosn tag '" + resource_options_tag + "'.");
	}
	resource_options_[ resource_options_tag ] = resource_options;
}


/////////////////////////////////////////////////////////
///// helper functions for resource configurations //////
/////////////////////////////////////////////////////////

ResourceConfiguration const &
LazyResourceManager::find_resource_configuration(
	ResourceTag const & resource_tag
) const {
	auto config(
		resource_configurations_.find(resource_tag));
	if ( config == resource_configurations_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "Unable to find resource configuration for the resource tag '" + resource_tag + "'.");
	}
	return config->second;
}

ResourceLocatorOP
LazyResourceManager::find_resource_locator(
	LocatorTag const & locator_tag
) const {
	auto resource_locator(
		resource_locators_.find(locator_tag));
	if ( resource_locator == resource_locators_.end() ) {
		throw utility::excn::EXCN_Msg_Exception("Unable to find resource locator for the locator tag '" + locator_tag + "'.");
	}
	return resource_locator->second;
}

ResourceOptionsOP
LazyResourceManager::find_resource_options(
	ResourceOptionsTag const & options_tag
) const {
	auto resource_options(
		resource_options_.find(options_tag));
	if ( resource_options == resource_options_.end() ) {
		throw utility::excn::EXCN_Msg_Exception("Unable to find resource options for the resource options tag '" + options_tag + "'.");
	}
	return resource_options->second;
}

bool
LazyResourceManager::has_resource(
	ResourceTag const & resource_tag
) const {
	if ( ResourceManager::has_resource(resource_tag) ) {
		return true;
	} else if ( has_resource_configuration(resource_tag) ) {
		return true;
	} else {
		return false;
	}
}

ResourceOP
LazyResourceManager::find_resource(
	ResourceTag const & resource_tag
) {
	if ( ResourceManager::has_resource(resource_tag) ) {
		return ResourceManager::find_resource(resource_tag);
	} else if ( has_resource_configuration(resource_tag) ) {
		ResourceOP new_resource(create_resource(resource_tag));
		add_resource(resource_tag, new_resource);
		return new_resource;
	} else {
		stringstream err_msg;
		err_msg
			<< "Unable to find resource with tag '" << resource_tag << "'." << endl;
		utility_exit_with_message(err_msg.str());
	}
}


ResourceOP
LazyResourceManager::create_resource(
	LocatorTag const & resource_tag
) const {

	ResourceConfiguration const & config(
		find_resource_configuration(resource_tag));

	ResourceLocatorOP locator(
		find_resource_locator(config.locator_tag));

	ResourceStreamOP stream(
		locator->locate_resource_stream(config.locator_id));

	ResourceLoaderOP loader(
		ResourceLoaderFactory::get_instance()->create_resource_loader(
		config.loader_type));

	ResourceOptionsOP resource_options;
	if ( config.resource_options_tag != "" ) {
		resource_options = find_resource_options(config.resource_options_tag);
	} else {
		resource_options = loader->default_options();
	}

	ResourceOP resource(
		loader->create_resource(
		*resource_options,
		config.locator_id,
		stream->stream()));

	return resource;
}

void
LazyResourceManager::create_resources(
	JobTag const &
) {
	utility_exit_with_message("This is meant to be overwritten in the the derived class.");
}

void
LazyResourceManager::show(
	std::ostream & out
) const {
	ResourceManager::show(out);
	out << endl;

	out
		<< "LazyResourceManager.resource_tags:" << endl
		<< std::setiosflags(std::ios::left) << setw(16) << "ResourceDescription"
		<< std::setiosflags(std::ios::left) << setw(16) << "JobTag ->"
		<< "ResourceTag" << endl;
	 for ( auto const & resource_tag : resource_tags_ ) {
		out
			<< std::setiosflags(std::ios::left) << setw(16) << resource_tag.first.first << setw(16) << resource_tag.first.second
			<< resource_tag.second << endl;
	}
	out << endl;

	out
		<< "LazyResourceManager.job_options:" << endl;
	 for ( auto const & job_option : job_options_ ) {
		out
			<< "ResourceTag: " << job_option.first << endl
			<< *(job_option.second) << endl;
	}
	out << endl;

	out
		<< "LazyResourceManager.resource_configurations:" << endl;
	 for ( auto const & resource_configuration : resource_configurations_ ) {
		out
			<< "ResourceTag: " << resource_configuration.first << endl
			<< resource_configuration.second << endl;
	}
	out << endl;

	out
		<< "LazyResourceManager.resource_locators:" << endl;
	 for ( auto const & resource_locator : resource_locators_ ) {
		out
			<< "LocatorTag: " << resource_locator.first << endl;
		resource_locator.second->show(out);
		out << endl;
	}
	out << endl;

	out
		<< "LazyResourceManager.resource_options:" << endl;
	 for ( auto const & resource_option : resource_options_ ) {
		out
			<< "LocatorOptionsTag: " << resource_option.first << endl
			<< *(resource_option.second) << endl;
	}

}


bool
LazyResourceManager::has_resource_configuration( ResourceTag const & resource_tag ) const
{
	return resource_configurations_.find( resource_tag ) != resource_configurations_.end();
}

bool
LazyResourceManager::has_resource_locator( LocatorTag const & locator_tag ) const
{
	return resource_locators_.find( locator_tag ) != resource_locators_.end();
}
bool
LazyResourceManager::has_resource_options( ResourceOptionsTag const & resource_options_tag ) const
{
	return resource_options_.find( resource_options_tag ) != resource_options_.end();
}

std::ostream &
operator<<(
	std::ostream & out,
	LazyResourceManager const & lazy_resource_manager
) {
	lazy_resource_manager.show(out);
	return out;
}


} // namespace
} // namespace
