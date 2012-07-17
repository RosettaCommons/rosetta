// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/LazyResourceManager.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_LazyResourceManager_hh
#define INCLUDED_basic_resource_manager_LazyResourceManager_hh

// Unit Headers
#include <basic/resource_manager/LazyResourceManager.fwd.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/JobOptions.hh>

// Project Headers
#include <basic/resource_manager/types.hh>
#include <basic/resource_manager/ResourceLoader.fwd.hh>
#include <basic/resource_manager/ResourceLocator.fwd.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <utility/exit.hh>

namespace basic {
namespace resource_manager {

struct ResourceConfiguration {
	ResourceTag resource_tag;
	LocatorTag locator_tag;
	LocatorID locator_id;
	LoaderType loader_type;
	ResourceOptionsTag resource_options_tag;
};

class LazyResourceManager : public ResourceManager
{
protected:
	LazyResourceManager();

public:
	virtual ~LazyResourceManager();


public: // work with resources by ResourceDescription + JobTag
	virtual
	ResourceOP
	create_resource_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag) const;

	virtual
	void
	add_resource_tag_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag,
		ResourceTag const & resource_tag);

	virtual
	bool
	has_resource_tag_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag) const;

	virtual
	ResourceTag
	find_resource_tag_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag) const;

	virtual
	ResourceOP
	get_resource_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag);

	virtual
	void
	free_resource_by_job_tag(
		ResourceDescription const & resource_description,
		JobTag const & job_tag);

	/// @brief wipe out all data; useful for unit testing, but hard to fathom how it would be useful
	/// otherwise.
	virtual
	void
	clear();

public: // Work with options with OptionKey + JobTag

	virtual
	void
	add_job_options(
		JobTag const & job_tag,
		JobOptionsOP job_options);

	virtual
	bool
	has_job_options(
		JobTag const & job_tag) const;

	virtual
	JobOptionsOP
	get_job_options(
		JobTag const & job_tag) const;


public: // Set resource configuration for each resource
	virtual
	void
	add_resource_configuration(
		ResourceTag const & resource_tag,
	  ResourceConfiguration const & resource_configuration);

	virtual
	void
	add_resource_locator(
		LocatorTag const & locator_tag,
		ResourceLocatorOP resource_locator);

	virtual
	void
	add_resource_options(
		ResourceOptionsTag const & resource_options_tag,
		ResourceOptionsOP resource_options);

public: // helper functions relating to resource configuration and creation
	virtual
	ResourceConfiguration const &
	find_resource_configuration(
		ResourceTag const & resource_tag) const;

	virtual
	ResourceLocatorOP
	find_resource_locator(
		LocatorTag const & locator_tag) const;

	virtual
	ResourceOptionsOP
	find_resource_options(
		ResourceOptionsTag const & resource_options_tag) const;

	virtual
	ResourceOP
	create_resource(
		ResourceTag const & resource_tag) const;

public: // Interface to for creating and accessing resources

	/// @brief Create all the resources for a particular job; this should be implemented by the derived
	/// class.  The implementation in this class calls utility::exit
	virtual
	void
	create_resources( JobTag const & );

	///@brief called by the protocol, return initialized resource; this should be implemented by the
	/// derived class.  The implementation in this class calls utility::exit
	virtual
	ResourceOP
	get_resource( ResourceDescription const & );

	/// @brief has a ResourceConfiguration been provided to the LazyResourceManager for a Resource with a particular ResourceTag?
	bool has_resource_configuration( ResourceTag const & resource_tag ) const;
	/// @brief has a ResourceLocator object been provided to the LazyResourceManager which has a particular LocatorTag?
	bool has_resource_locator( LocatorTag const & locator_tag ) const;
	/// @brief has a ResourceOptions object been provided to the LazyResourceManager which has a particular ResourceOptionsTag?
	bool has_resource_options( ResourceOptionsTag const & resource_options_tag ) const;

private: // Data members

	typedef std::map< std::pair< ResourceDescription, JobTag >,	ResourceTag >
		ResourceTagsMap;
	ResourceTagsMap resource_tags_;

	typedef std::map< ResourceTag, JobOptionsOP >
		JobOptionsMap;
	JobOptionsMap job_options_;

	typedef std::map< ResourceTag, ResourceConfiguration >
		ResourceConfigurationMap;
	ResourceConfigurationMap resource_configurations_;

	typedef std::map< LocatorTag, ResourceLocatorOP > ResourceLocatorsMap;
	ResourceLocatorsMap resource_locators_;

	typedef std::map< ResourceOptionsTag, ResourceOptionsOP >
		ResourceOptionsMap;
	ResourceOptionsMap resource_options_;

};

} // namespace resource_manager
} // namespace basic

#endif
