// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManager.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManager_hh
#define INCLUDED_basic_resource_manager_ResourceManager_hh

//unit headers
#include <basic/resource_manager/ResourceManager.fwd.hh>
#include <basic/resource_manager/types.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.fwd.hh>
#include <basic/resource_manager/ResourceLocator.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/XMLSchemaValidation.fwd.hh>

//C++ headers
#include <list>
#include <map>
#include <set>
#include <string>

namespace basic {
namespace resource_manager {

class ResourceManager : public utility::pointer::ReferenceCount
{
public:

	ResourceManager();
	~ResourceManager() override;

	/// @brief Read all of the resource definition files specified on the command line
	/// from the "-resource_definitions <fname1> <fname2>... " flag.
	void
	initialize_from_commandline();

	/// @brief The ResourceManager will take the input stream and
	/// parse it as an XML document. Then it will read the ResourceLocators
	/// and Resouces tags out of this stream. Multiple streams may be read
	/// this way, though, if two streams try and declare the locators or resources
	/// that have the same name, then the ResourceManager will throw an exception
	void read_resources_from_xml( std::string const & filename, std::string const & input_xml );

	/// @brief How many resource definitions were read into the %ResourceManager?
	platform::Size
	n_resources_declared() const;

	/// @brief How many resources are currently being held by the %ResourceManager?
	platform::Size
	n_resources_in_memory() const;

	/// @breif Has a particular resource been declared?
	bool
	has_resource(
		std::string const & resource_name
	) const;

	/// @breif Does a particular resource currently reside in memory? Mostly this
	/// function is for testing purposes.
	bool
	has_resource_in_memory(
		std::string const & resource_name
	) const;

	/// @brief Return a list of all of the resources that have been declared.
	std::list< std::string >
	resources_that_have_been_declared() const;

	/// @brief When a resource no longer needs to be held in memory, it may be
	/// deallocated. If the resource has not been allocated, this function
	/// simply returns (without throwing an exception).
	void
	deallocate_resource(
		std::string const & resource_name );

	/// @brief Delete all stored data, including all information about how to construct
	/// not-yet-constructed resources.
	void
	clear();

	/// @brief Get a resource with a given name.
	ResourceCOP
	get_resource(
		std::string const & resource_name );

	/// @brief Construct the XML Schema for all of the various resource types that can be
	/// loaded using the ResourceManager using all of the various input-stream-fetching
	/// techniques (i.e. file on disk vs database query).
	static
	std::string
	schema_for_resource_definition_file();

	/// @brief The name mangling function for the top-level elements of the <ResourceDefinitions>
	/// tag in the resource-definition file.
	static
	std::string
	complex_type_name_for_resource_def( std::string const & element_name );

private:

	/// @brief After having loaded a given XML file in from disk, validate its contents against
	/// the internally generated XML schema.
	void
	validate_input_against_xsd( std::string const & xml_filename, std::string const & input_xml );

	/// @brief Process all of the subtags of the given <Resources> tag
	void
	read_resources_tags(
		std::string const & filename,
		utility::tag::TagCOP tags );

	/// @brief Read the portion of an XML file that declares ResourceLocator objects
	void
	read_resource_locators_tags(
		std::string const & filename,
		utility::tag::TagCOP tags );

	/// @brief Create and store a particular Resource that has been requested and that has not yet
	/// been constructed.
	ResourceCOP
	create_resource(
		std::string const & resource_name );

	/// @brief Create and store a particular ResourceLocator that has been requested and that
	/// has not yet been constructed
	ResourceLocatorOP
	create_resource_locator(
		std::string const & resource_locator_name ) const;

	/// @brief Return the ResourceLocator for the given resource, potentially instantiating
	/// this locator if it has not been previously instantiated.
	ResourceLocatorOP
	find_resource_locator(
		std::string const & locator_name );

	std::string
	locator_name_from_resource_tag(
		utility::tag::TagCOP resource_tag
	) const;

	std::string
	input_id_from_resource_tag(
		utility::tag::TagCOP resource_tag
	) const;

	std::string
	loader_type_from_resource_tag(
		utility::tag::TagCOP resource_tag
	) const;

	/// @brief Do not deallocate the indicated resource until explicitly requested.
	/// This is only necessary for resources that are loaded during the course of
	/// constructing another resource, and that would be deallocated automatically
	/// when the first resource was deallocated.
	void
	require_explicit_deallocation_for_resource( std::string const & resource_name );

private:

	typedef std::map< std::string, std::string > DeclarationMap;
	/// @brief The file in which a particular resource locator was declared; useful if two resource locators
	/// from two different files were given the same name.
	DeclarationMap resource_locator_declarations_;
	/// @brief The file in which a particular resource was declared; useful if two resources
	/// from two different files were given the same name.
	DeclarationMap resource_declarations_;


	/// @brief The instructions on how to create each resource locator that has been declared
	typedef std::map< std::string, utility::tag::TagCOP > ResourceLocatorTagMap;
	ResourceLocatorTagMap resource_locator_tags_;

	/// @brief The instructions on how to create each resource that has been declared
	typedef std::map< std::string, utility::tag::TagCOP > ResourceTagMap;
	ResourceTagMap resource_tags_;

	/// @brief The locator to use in order to build a particular resource
	/// If a resource does not have an entry in this map, then the default
	/// FileSysteResoureLocator should be used.
	typedef std::map< std::string, std::string > LocatorsForResourcesMap;
	LocatorsForResourcesMap locators_for_resources_;

	/// @brief The set of Locators
	typedef std::map< std::string, ResourceLocatorOP > LocatorMap;
	LocatorMap resource_locators_;

	typedef std::map< std::string, ResourceCOP > ResourcesMap;
	ResourcesMap resources_;

	/// @brief The object used to validate XML input files against the schema
	utility::tag::XMLSchemaValidatorOP validator_;

	////// Data to track the interdependencies of resources on each other: one resource may request another
	////// Resource as it is being constructed. When this happens, the resource manager records the dependency
	////// between these two resources so that later, when the first resource is deallocated, it can also
	////// deallocate the second -- but only so long as no other resource also depends on the second one.
	utility::vector1< std::string > resources_being_constructed_stack_;
	std::map< std::string, std::set< std::string > > resources_a_resource_depends_on_;
	std::map< std::string, std::set< std::string > > living_resources_dependent_on_a_resource_;
	std::set< std::string > resources_to_hold_until_deallocation_explicitly_requested_;
};

} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_hh
