// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManager.cc
/// @brief
/// @author


//unit headers
#include <basic/resource_manager/ResourceManager.hh>

//project headres
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/ResourceLocatorCreator.hh>
#include <basic/resource_manager/ResourceLocatorFactory.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>

//utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/thread/threadsafe_creation.hh>

//basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

//C++ headers
#include <sstream>
#include <iomanip>

namespace basic {
namespace resource_manager {

//using std::stringstream;
//using std::endl;
//using std::setw;

static Tracer TR( "basic.resource_manager.ResourceManager" );

ResourceManager::ResourceManager() :
	validator_( new utility::tag::XMLValidator )
{}

ResourceManager::~ResourceManager() = default;

void
ResourceManager::initialize_from_commandline()
{
	utility::vector1< std::string > resdef_files = basic::options::option[ basic::options::OptionKeys::jd3::resource_definition_files ]();
	for ( std::string const & fname : resdef_files ) {
		utility::io::izstream ifs( fname.c_str() );
		if ( ! ifs.good() ) {
			std::ostringstream oss;
			oss << "Error in trying to read the resources definition file '" << fname << "'; file could not be found.\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
		std::string file_contents; utility::slurp( ifs, file_contents );
		read_resources_from_xml( fname, file_contents );
	}
}

void
ResourceManager::read_resources_from_xml(
	std::string const & filename,
	std::string const & input_xml
)
{
	validate_input_against_xsd( filename, input_xml );
	utility::tag::TagCOP tags = utility::tag::Tag::create( input_xml );
	for ( auto const & subtag : tags->getTags() ) {
		debug_assert( subtag->getName() == "Resources" || subtag->getName() == "ResourceLocators" );
		if ( subtag->getName() == "Resources" ) {
			read_resources_tags( filename, subtag );
		} else {
			read_resource_locators_tags( filename, subtag );
		}
	}
}

/// @brief How many resource definitions were read into the %ResourceManager?
platform::Size
ResourceManager::n_resources_declared() const
{
	return resource_tags_.size();
}

/// @brief How many resources are currently being held by the %ResourceManager?
platform::Size
ResourceManager::n_resources_in_memory() const
{
	return resources_.size();
}


bool
ResourceManager::has_resource(
	std::string const & resource_name
) const
{
	return resource_tags_.count( resource_name );
}


bool
ResourceManager::has_resource_in_memory(
	std::string const & resource_name
) const
{
	return resources_.count( resource_name ) != 0;
}


std::list< std::string >
ResourceManager::resources_that_have_been_declared() const
{
	std::list< std::string > resource_names;
	for ( auto const & resource_pair : resource_tags_ ) {
		resource_names.push_back( resource_pair.first );
	}
	return resource_names;
}

void
ResourceManager::deallocate_resource(
	std::string const & resource_name
) {
	auto iter = resources_.find( resource_name );
	if ( iter != resources_.end() ) {
		resources_.erase( iter );
	}
	if ( resources_to_hold_until_deallocation_explicitly_requested_.count( resource_name ) ) {
		// if a resource is being explicitly deleted, then, if another resource comes along
		// and requests it to be re-instantiated, it should not remain in memory indefinitely.
		// remove it from this set.
		resources_to_hold_until_deallocation_explicitly_requested_.erase( resource_name );
	}

	// recursively deallocate the resources that are dependent on this resource
	// and that have not been requested to remain in memory until they are explicitly
	// deleted
	if ( resources_a_resource_depends_on_.count( resource_name ) ) {
		for ( auto const & dep_res : resources_a_resource_depends_on_[ resource_name ] ) {
			debug_assert( living_resources_dependent_on_a_resource_.count( dep_res ) );
			debug_assert( living_resources_dependent_on_a_resource_[ dep_res ].count( resource_name ) );
			living_resources_dependent_on_a_resource_[ dep_res ].erase( resource_name );
			if ( living_resources_dependent_on_a_resource_[ dep_res ].empty() &&
					! resources_to_hold_until_deallocation_explicitly_requested_.count( dep_res ) ) {
				living_resources_dependent_on_a_resource_.erase( dep_res );
				deallocate_resource( dep_res );
			}
		}
		resources_a_resource_depends_on_.erase( resource_name );
	}
}

void
ResourceManager::clear()
{
	resource_locator_declarations_.clear();
	resource_declarations_.clear();
	resource_locator_tags_.clear();
	resource_tags_.clear();
	locators_for_resources_.clear();
	resource_locators_.clear();
	resources_.clear();
	resources_a_resource_depends_on_.clear();
	living_resources_dependent_on_a_resource_.clear();
	resources_to_hold_until_deallocation_explicitly_requested_.clear();
}

ResourceCOP
ResourceManager::get_resource(
	std::string const & resource_name
)
{
	for ( std::string const & being_constructed_resource_name : resources_being_constructed_stack_ ) {
		if ( being_constructed_resource_name == resource_name ) {
			std::ostringstream oss;
			oss << "Critical error: a resource may not depend on another resource that"
				" in turn depends on the first one. Infinite recursion is not allowed.\n";
			oss << "Resources that were being created were:\n";
			for ( std::string const & rn : resources_being_constructed_stack_ ) {
				oss << "   " << rn << "\n";
			}
			oss << "    ultimately requesting resource " << resource_name << "\n";
			oss << "Cannot proceed from here.\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
	}

	if ( resources_being_constructed_stack_.size() != 0 ) {
		resources_a_resource_depends_on_[ resources_being_constructed_stack_.back() ].insert( resource_name );
		living_resources_dependent_on_a_resource_[ resource_name ].insert( resources_being_constructed_stack_.back() );
	}

	auto iter = resources_.find( resource_name );
	if ( iter == resources_.end() ) {
		if ( ! has_resource( resource_name ) ) {
			std::ostringstream oss;
			oss << "Cannot create the requested resource: '" << resource_name << "'; has this resource been declared?\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
		return create_resource( resource_name );
	}
	return iter->second;
}

std::string
ResourceManager::schema_for_resource_definition_file()
{
	TR << "Generating schema for ResourceManager..." << std::endl;
	using namespace utility::tag;
	XMLSchemaDefinition xsd;

	// Example: (Ignore the period on the left; it's to prevent re-indentation by the beautifier)
	//. <ResourceDefinitions>
	//.   <ResourceLocators>
	//.     <FileSystemLocator name="std" search_paths="/home/andrew/top4400/:/home/andrew/msd/inputs/" />
	//.     <FileSystemLocator name="designs_dir" search_paths="/home/andrew/designs/des001/" />
	//.     <DatabaseLocator name="db" file="features.db3" table="structures" />
	//.   </ResourceLocators>
	//.   <Resources>
	//.     <Resource name="1l2y_native" input_id="1l2yH_nowater.pdb" locator="std">
	//.       <PoseResource>
	//.         <PDB>
	//.           <PoseLoaderOptions skip_conect_info="false"/>
	//.         <PDB>
	//.       </PoseResource>
	//.     </Resource>
	//.     <Resource name="2ten_homologue_1_cen" input_id="2ten_binary.out">
	//.       <PoseResource>
	//.         <SilentFile silent_tag="2ten_0001">
	//.           <SilentFileOptions in_fullatom="false"/>
	//.         </SilentFile>
	//.       </PoseResource>
	//.     </Resource>
	//.     <Resource name="1l2y_3mers_exclude_homologues_2014db" input_id="1l2y_nohom14.frags3">
	//.       <FragmentSet ...options for a fragment set?... >
	//.     </Resource>
	//.   </Resources>
	//. </ResourceDefiitions>

	std::string const resources_name = "Resources";
	std::string const resource_name = "Resource";
	//std::string const resource_loaders_name = "ResourceLoaders";
	std::string const resource_locators_name = "ResourceLocators";
	std::string const resource_definitions_name = "ResourceDefinitions";

	XMLSchemaElement resource_definitions_element;
	resource_definitions_element.name( resource_definitions_name )
		.type_name( complex_type_name_for_resource_def( resource_definitions_name ));
	xsd.add_top_level_element( resource_definitions_element );

	// ResourceDefinitions complex type:
	XMLSchemaSimpleSubelementList resource_definitions_subelements;
	resource_definitions_subelements.add_already_defined_subelement( resources_name, &complex_type_name_for_resource_def );
	resource_definitions_subelements.add_already_defined_subelement( resource_locators_name, &complex_type_name_for_resource_def );

	XMLSchemaComplexTypeGenerator resource_definitions_ct;
	resource_definitions_ct.element_name( resource_definitions_name )
		.complex_type_naming_func( &complex_type_name_for_resource_def )
		.description( "The ResourceDefinitions block contains as many ResourceLocators and Resources sub-tags as desired"
		" in any order. A Resource is a large piece of data that is perfectly (bitwise) constant after its creation."
		" Resources can be initialized and they can be deallocated, but they cannot be modified. Because they are constant,"
		" Resources can be shared between two threads without worry about race condidtions. ResourceLocators are"
		" objects that help obtain the data stream that the Resource should be constructed from -- the default ResourceLocator"
		" is the FileSystem locator and it simply opens a file from the file system. Other more complex Locators can be imagined"
		" such as a DatabaseLocator which reads its data stream from some database perhaps over a network connection."
		" The purpose of the ResourceManager is to manage the loading and unloading of this data, which perhaps is expensive"
		" and should only be performed once."
		" Multiple resource definition files may be used to declare all of the resources that are provided;"
		" however, you cannot declare two resources or two resource locators with the same name in the set of files that"
		" you provide. If you have two resources that depend on the same locator, and it does not make sense to put the"
		" resources in the same file (e.g. perhaps it makes sense to have all of your PoseResources in one file, and all"
		" of your FragmentSetResources in another file) it is possible to create a third file with the common locator."
		" That is, it is possible to declare a Resource that depends on a ResourceLocator declared in a different file."
		" If a Resource requires a particular locator, that locator must be declared either earlier in the file,"
		" or in a file that is parsed before the file that the Resource is." )
		.set_subelements_repeatable( resource_definitions_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList resource_locators_subelements;
	for ( auto const & resource_locator_pair : ResourceLocatorFactory::get_instance()->locator_map() ) {
		resource_locators_subelements.add_already_defined_subelement(
			resource_locator_pair.first, &ResourceLocatorFactory::complex_type_name_for_locator );
		resource_locator_pair.second->provide_xml_schema( xsd );
	}

	XMLSchemaComplexTypeGenerator resource_locators_ct;
	resource_locators_ct.element_name( resource_locators_name )
		.complex_type_naming_func( & complex_type_name_for_resource_def )
		.description( "This block declares the set of ResourceLocators that are used to open streams for resource creation" )
		.set_subelements_repeatable( resource_locators_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList resources_subelements;
	resources_subelements.add_already_defined_subelement( resource_name, & complex_type_name_for_resource_def );

	XMLSchemaComplexTypeGenerator resources_ct;
	resources_ct.element_name( resources_name )
		.complex_type_naming_func( & complex_type_name_for_resource_def )
		.description( "In this block, Resources may be delcared" )
		.set_subelements_repeatable( resources_subelements )
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList resource_subelements;
	for ( auto const & resource_loader_pair : ResourceLoaderFactory::get_instance()->loader_map() ) {
		resource_subelements.add_already_defined_subelement(
			resource_loader_pair.first, &ResourceLoaderFactory::complex_type_name_for_loader );
		resource_loader_pair.second->provide_xml_schema( xsd );
	}

	typedef XMLSchemaAttribute Attr;
	AttributeList resource_attributes;
	resource_attributes
		+ required_name_attribute( "The name that will be used for this resource; must be unique among all resources that the ResourceManager reads" )
		+ Attr::required_attribute( "input_id", xs_string, "The string used by the ResourceLocator to obtain the stream that the resource should be read from; e.g. a file name" )
		+ Attr( "locator", xs_string, "The name of the ResourceLocator that should be used to load this resource" );

	XMLSchemaComplexTypeGenerator resource_ctgen;
	resource_ctgen.element_name( resource_name )
		.complex_type_naming_func( & complex_type_name_for_resource_def )
		.description( "The Resource tag declares a single resource. The Resource must be given a (unique) name and its"
		" input_id must be given so that the ResourceLocator can locate the right stream to open. If needed, a locator"
		" may be given, but a default locator -- the FileSystemLocator -- will be used if no locator is given. The"
		" Resource tag must have exactly one subtag which specifies the kind of Resource that is being created." )
		.set_subelements_pick_one( resource_subelements )
		.add_attributes( resource_attributes )
		.write_complex_type_to_schema( xsd );

	return xsd.full_definition();
}


std::string
ResourceManager::complex_type_name_for_resource_def( std::string const & element_name )
{
	return "resource_manager_" + element_name + "_type";
}

void
ResourceManager::validate_input_against_xsd( std::string const & xml_filename, std::string const & input_xml )
{
	// Step 1: Internal generation of the schema and initialization of the libxml2-schema-validation machinery
	if ( ! validator_->schema_has_been_set() ) {
		std::ostringstream oss;
		oss << "If you are seeing this message, the interanlly-generated XML Schema for the ResourceManager could not" <<
			" be properly generated\nThis failure occurred before the XML resource-definition file that was provided was" <<
			" examined. The error has been compiled into Rosetta and will need to be fixed by a developer.\n";
		std::string schema;
		try {
			schema = schema_for_resource_definition_file();
		} catch ( utility::excn::Exception const & e ) {
			oss << "An error was encountered while the string of the schema was being generated; this occurred before the schema was analyzed for whether it was correct.\n";
			oss << e.msg() << "\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}

		utility::tag::XMLValidationOutput schema_valid_output;
		try {
			TR << "Initializing schema validator..." << std::endl;
			schema_valid_output = validator_->set_schema( schema );
			TR << "...done" << std::endl;
		} catch ( utility::excn::Exception const & e ) {
			oss << e.msg() << "\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}

		if ( ! schema_valid_output.valid() ) {
			oss << "If there is an error message immediately above stating that the schema failed to validate, and errors below that like like 'real XML' with lots of <xs:something> tags and NOT like your XML input, then you have a global schema validation error and not an XML input validation error. Read the error message below and fix your schema in the C++ code.\n\n";
			oss << "Errors: " << schema_valid_output.error_messages() << "\n";
			oss << "Warnings: " << schema_valid_output.warning_messages() << "\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}
	}

	TR << "Validating input resource definition file..." << std::endl;
	std::ostringstream oss;
	oss << "Error in trying to validate the XML from the file '" << xml_filename << ". Input XML does not match the" <<
		" schema for a resource definition file. Error messages from the schema validator below:\n";


	utility::tag::XMLValidationOutput validator_output;
	try {
		validator_output = validator_->validate_xml_against_schema( input_xml );
	} catch ( utility::excn::Exception const & e ) {
		oss << e.msg() << "\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	if ( ! validator_output.valid() ) {
		oss << "Error messages were:\n" << validator_output.error_messages();
		oss << "------------------------------------------------\n";
		oss << "Warning messages were:\n" << validator_output.warning_messages();
		oss << "------------------------------------------------\n";
		throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	TR << "...done" << std::endl;
}

void
ResourceManager::read_resources_tags(
	std::string const & filename,
	utility::tag::TagCOP tag
)
{
	for ( auto const & subtag : tag->getTags() ) {
		debug_assert( subtag->hasOption( "name" ) );
		std::string resource_name = subtag->getOption< std::string >( "name" );
		if ( resource_tags_.count( resource_name ) ) {
			std::ostringstream oss;
			oss << "Error in reading the declaration of resource '" << resource_name << "' in the file '" << filename <<
				"' as a resource with the same name was declared previously in file '" <<
				resource_declarations_ << "'.\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}

		if ( subtag->hasOption( "locator" ) ) {
			std::string locator_name = subtag->getOption< std::string >( "locator" );

			if ( ! resource_locator_tags_.count( locator_name ) ) {
				std::ostringstream oss;
				oss << "Error in reading the declaration of resource '" << resource_name << "' in the file '" << filename <<
					"'. The ResourceLocator it has requested with the name '" << locator_name << "' has not been declared.\n";
				throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
			}
		}
		resource_tags_[ resource_name ] = subtag;
		resource_declarations_[ resource_name ] = filename;
	}
}

/// @brief Read the portion of an XML file that declares ResourceLocator objects
void
ResourceManager::read_resource_locators_tags(
	std::string const & filename,
	utility::tag::TagCOP resource_locators_tag
)
{
	for ( auto const & subtag : resource_locators_tag->getTags() ) {
		std::string locator_name = subtag->getOption< std::string >( "name" );
		if ( resource_locator_tags_.count( locator_name ) ) {
			std::ostringstream oss;
			oss << "Error in reading the ResourceLocator named '" << locator_name << "' in file '" << filename;
			oss << "' as this resource locator was previously declared in file '";
			oss << resource_locator_declarations_[ locator_name ] << "'.\n";
			throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
		}

		resource_locator_declarations_[ locator_name ] = filename;
		resource_locator_tags_[ locator_name ] = subtag;
	}
}

/// @brief Create and store a particular Resource that has been requested and that has not yet
/// been constructed.
ResourceCOP
ResourceManager::create_resource(
	std::string const & resource_name
)
{
	resources_being_constructed_stack_.push_back( resource_name );

	auto resource_tag_iter = resource_tags_.find( resource_name );
	debug_assert( resource_tag_iter != resource_tags_.end() );
	utility::tag::TagCOP resource_tag( resource_tag_iter->second );

	std::string locator_name = locator_name_from_resource_tag( resource_tag );
	ResourceLocatorOP locator( find_resource_locator( locator_name ) );

	std::string input_id = input_id_from_resource_tag( resource_tag);
	ResourceStreamOP stream( locator->locate_resource_stream( input_id ));

	ResourceLoaderOP loader(
		ResourceLoaderFactory::get_instance()->create_resource_loader(
		loader_type_from_resource_tag( resource_tag ) ));

	ResourceCOP resource( loader->create_resource(
		*this,
		resource_tag->getTags()[0],
		input_id,
		stream->stream()));

	resources_[ resource_name ] = resource;
	resources_being_constructed_stack_.pop_back();
	return resource;
}

/// @brief Create and store a particular ResourceLocator that has been requested and that
/// has not yet been constructed
ResourceLocatorOP
ResourceManager::create_resource_locator(
	std::string const & resource_locator_name
) const
{
	ResourceLocatorOP locator;

	auto tag_iter = resource_locator_tags_.find( resource_locator_name );
	if ( tag_iter != resource_locator_tags_.end() ) {
		locator = ResourceLocatorFactory::get_instance()->create_resource_locator(
			tag_iter->second->getName(), resource_locator_name, tag_iter->second );
	} else if ( resource_locator_name == "" ) {
		locator.reset( new locator::FileSystemResourceLocator );
	} else {
		// ERROR!!!
	}
	return locator;
}

/// @brief Return the ResourceLocator for the given resource, potentially instantiating
/// this locator if it has not been previously instantiated.
ResourceLocatorOP
ResourceManager::find_resource_locator(
	std::string const & locator_name
)
{
	auto locator_iter = resource_locators_.find( locator_name );
	if ( locator_iter != resource_locators_.end() ) {
		return locator_iter->second;
	}

	auto locator = create_resource_locator( locator_name );
	resource_locators_[ locator_name ] = locator;
	return locator;
}

std::string
ResourceManager::locator_name_from_resource_tag(
	utility::tag::TagCOP resource_tag
) const
{
	if ( resource_tag->hasOption( "locator" ) ) {
		return resource_tag->getOption< std::string >( "locator" );
	} else {
		// If no locator is given, then the default locator -- the file system locator --
		// is implicitly being requested. It is designated by the empty string.
		return "";
	}
}

std::string
ResourceManager::input_id_from_resource_tag(
	utility::tag::TagCOP resource_tag
) const
{
	return resource_tag->getOption< std::string >( "input_id" );
}

std::string
ResourceManager::loader_type_from_resource_tag(
	utility::tag::TagCOP resource_tag
) const
{
	debug_assert( resource_tag->getTags().size() > 0 );
	utility::tag::TagCOP subtag = resource_tag->getTags()[ 0 ];
	return subtag->getName();
}

void
ResourceManager::require_explicit_deallocation_for_resource( std::string const & resource_name )
{
	resources_to_hold_until_deallocation_explicitly_requested_.insert( resource_name );
}


} // namespace resource_manager
} // namespace basic
