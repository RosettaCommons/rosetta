// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/resource_manager/ResourceManager.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceManager
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <basic/resource_manager/ResourceLoaderRegistrator.hh>
#include <basic/resource_manager/loader_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>

//#include <basic/resource_manager/;

using namespace basic::resource_manager;

class DummyResource2 : public utility::pointer::ReferenceCount
{
public:
	DummyResource2( int value ) : value_( value ) {}
	int value_;
};

class DummyResource2Loader : public basic::resource_manager::ResourceLoader
{
public:
	DummyResource2Loader() {}

	static
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
		using namespace utility::tag;
		AttributeList attlist;
		attlist + XMLSchemaAttribute( "subresources", xs_string, "request another resource when constructing this one" );
		basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd, "Dummy", "(desc)", attlist );
	}

	basic::resource_manager::ResourceCOP
	create_resource(
		basic::resource_manager::ResourceManager & rm,
		utility::tag::TagCOP tag,
		std::string const &,
		std::istream &
	) const override
	{
		if ( tag->hasOption( "subresources" ) ) {
			// perhaps will throw!
			utility::vector1< std::string > subresources = utility::string_split( tag->getOption<std::string>( "subresources" ), ',' );
			for ( std::string const & subres : subresources ) {
				rm.get_resource( subres );
			}
		}
		return std::make_shared< DummyResource2 >( ++count_resources );
	}

	static int count_resources;

};

int
DummyResource2Loader::count_resources( 0 );


class DummyResource2LoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	/// @brief Return a up-casted owning pointer (ResourceLoaderOP) to the resource loader.
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const override {
		return std::make_shared< DummyResource2Loader >();
	}

	/// @brief Return the string identifier for the associated ResourceLoader
	std::string loader_type() const override {
		return "Dummy";
	}

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override
	{
		DummyResource2Loader::provide_xml_schema( xsd );
	}

};

class ResourceManagerTests : public CxxTest::TestSuite {
private:
	bool registered_dummy_ = { false };

public:

	void setUp() {
		if ( ! registered_dummy_ ) {
			registered_dummy_ = true;
			basic::resource_manager::ResourceLoaderRegistrator< DummyResource2LoaderCreator > reg;
		}
		core_init();
	}

	void test_schema_structure() {
		TS_ASSERT( true );
		std::string schema = ResourceManager::schema_for_resource_definition_file();
		utility::tag::TagCOP schema_tag( utility::tag::Tag::create( schema ) );

		TS_ASSERT( schema_tag->hasTag( "xs:element" ) );
		TS_ASSERT( schema_tag->hasTag( "xs:complexType" ) );
		utility::tag::TagCOP resource_definitions_element_tag;
		utility::tag::TagCOP resource_definitions_ct_tag;
		std::string resource_definitions_ct_name( ResourceManager::complex_type_name_for_resource_def( "ResourceDefinitions" ) );

		for ( auto const & subtag : schema_tag->getTags( "xs:element" ) ) {
			// there really only should be a single element subtag
			if ( subtag->hasOption( "name" ) && subtag->getOption< std::string >( "name" ) == "ResourceDefinitions" ) {
				resource_definitions_element_tag = subtag;
				break;
			}
		}

		for ( auto const & subtag : schema_tag->getTags( "xs:complexType" ) ) {
			std::string subtag_name_attr = subtag->getOption< std::string >( "name", "" );
			if ( subtag_name_attr == resource_definitions_ct_name ) {
				resource_definitions_ct_tag = subtag;
				break;
			}
		}


		TS_ASSERT( resource_definitions_element_tag );
		{
			// there should be two subelements of the ResourceDefinitions element
			utility::tag::TagCOP rde_ct_tag = resource_definitions_ct_tag->getTag( "xs:choice" );
			TS_ASSERT( rde_ct_tag );
			bool resources_found( false ), resource_locators_found( false );
			for ( auto const & subtag : rde_ct_tag->getTags( "xs:element" ) ) {
				if ( subtag->getOption< std::string >( "name", "" ) == "Resources" ) {
					resources_found = true;
				} else if ( subtag->getOption< std::string >( "name", "" ) == "ResourceLocators" ) {
					resource_locators_found = true;
				}
			}
			TS_ASSERT( resources_found );
			TS_ASSERT( resource_locators_found );
		}

		{
			// Now let's look at the subtags of the top-level tag and see if we can't find the
			// complex types for the Resources and ResourceLocators tags
			std::string resources_ct_name( ResourceManager::complex_type_name_for_resource_def( "Resources" ) );
			std::string resource_locators_ct_name( ResourceManager::complex_type_name_for_resource_def( "ResourceLocators" ) );
			bool found_resources_ct( false ), found_resource_locators_ct( false );
			for ( auto const & subtag : schema_tag->getTags( "xs:complexType" ) ) {
				std::string subtag_name_attr = subtag->getOption< std::string >( "name", "" );
				if ( subtag_name_attr == resources_ct_name ) {
					found_resources_ct = true;
				} else if ( subtag_name_attr == resource_locators_ct_name ) {
					found_resource_locators_ct = true;
				}
				if ( found_resources_ct && found_resource_locators_ct ) break;
			}
			TS_ASSERT( found_resources_ct );
			TS_ASSERT( found_resource_locators_ct );
		}

	}

	void test_parse_resource_definition_file() {
		using namespace basic::resource_manager;

		// Define two resources
		std::string resource_def_file =
			"<ResourceDefinitions>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator name=\"unit_test_files\" search_paths=\"protocols/rosetta_scripts/\"/>\n"
			"  </ResourceLocators>\n"
			"  <Resources>\n"
			"    <Resource name=\"1abc\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy/>\n"
			"    </Resource>\n"
			"    <Resource name=\"2def\" input_id=\"2def.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy/>\n"
			"    </Resource>\n"
			"  </Resources>\n"
			"</ResourceDefinitions>\n";

		ResourceManager rm;
		rm.read_resources_from_xml( "(unit test)", resource_def_file );
		TS_ASSERT_EQUALS( rm.n_resources_declared(), 2 );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		rm.get_resource( "1abc" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 1 );

		rm.deallocate_resource( "1abc" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		rm.get_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 1 );

		rm.deallocate_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );
	}

	void test_parse_resource_definition_w_dependency() {
		using namespace basic::resource_manager;

		// Define two resources
		std::string resource_def_file =
			"<ResourceDefinitions>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator name=\"unit_test_files\" search_paths=\"protocols/rosetta_scripts/\"/>\n"
			"  </ResourceLocators>\n"
			"  <Resources>\n"
			"    <Resource name=\"1abc\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"2def\"/>\n"
			"    </Resource>\n"
			"    <Resource name=\"2def\" input_id=\"2def.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy/>\n"
			"    </Resource>\n"
			"  </Resources>\n"
			"</ResourceDefinitions>\n";

		ResourceManager rm;
		rm.read_resources_from_xml( "(unit test)", resource_def_file );
		TS_ASSERT_EQUALS( rm.n_resources_declared(), 2 );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		// This should cause both 1abc and 2def to be loaded into memory
		rm.get_resource( "1abc" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 2 );

		rm.deallocate_resource( "1abc" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		rm.get_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 1 );

		rm.deallocate_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );
	}

	void test_parse_resource_definition_file_w_multiple_dependencies() {
		using namespace basic::resource_manager;

		// Define two resources
		std::string resource_def_file =
			"<ResourceDefinitions>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator name=\"unit_test_files\" search_paths=\"protocols/rosetta_scripts/\"/>\n"
			"  </ResourceLocators>\n"
			"  <Resources>\n"
			"    <Resource name=\"1abc\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy/>\n"
			"    </Resource>\n"
			"    <Resource name=\"2def\" input_id=\"2def.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"1abc\"/>\n"
			"    </Resource>\n"
			"    <Resource name=\"3ghi\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"1abc\"/>\n"
			"    </Resource>\n"
			"  </Resources>\n"
			"</ResourceDefinitions>\n";

		ResourceManager rm;
		rm.read_resources_from_xml( "(unit test)", resource_def_file );
		TS_ASSERT_EQUALS( rm.n_resources_declared(), 3 );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		rm.get_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 2 );
		TS_ASSERT( rm.has_resource_in_memory( "2def" ) );
		TS_ASSERT( rm.has_resource_in_memory( "1abc" ) );

		rm.get_resource( "3ghi" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 3 );
		TS_ASSERT( rm.has_resource_in_memory( "3ghi" ) );
		TS_ASSERT( rm.has_resource_in_memory( "2def" ) );
		TS_ASSERT( rm.has_resource_in_memory( "1abc" ) );

		rm.deallocate_resource( "2def" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 2 );
		TS_ASSERT( rm.has_resource_in_memory( "3ghi" )   );
		TS_ASSERT( ! rm.has_resource_in_memory( "2def" ) );
		TS_ASSERT( rm.has_resource_in_memory( "1abc" )   );

		rm.deallocate_resource( "3ghi" );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

	}

	void test_parse_resource_definition_file_w_circular_dependencies() {
		using namespace basic::resource_manager;

		// Define two resources
		std::string resource_def_file =
			"<ResourceDefinitions>\n"
			"  <ResourceLocators>\n"
			"    <FileSystemResourceLocator name=\"unit_test_files\" search_paths=\"protocols/rosetta_scripts/\"/>\n"
			"  </ResourceLocators>\n"
			"  <Resources>\n"
			"    <Resource name=\"1abc\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"2def\"/>\n"
			"    </Resource>\n"
			"    <Resource name=\"2def\" input_id=\"2def.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"3ghi\"/>\n"
			"    </Resource>\n"
			"    <Resource name=\"3ghi\" input_id=\"1abc.pdb\" locator=\"unit_test_files\">\n"
			"      <Dummy subresources=\"1abc\"/>\n"
			"    </Resource>\n"
			"  </Resources>\n"
			"</ResourceDefinitions>\n";

		ResourceManager rm;
		rm.read_resources_from_xml( "(unit test)", resource_def_file );
		TS_ASSERT_EQUALS( rm.n_resources_declared(), 3 );
		TS_ASSERT_EQUALS( rm.n_resources_in_memory(), 0 );

		try {
			rm.get_resource( "1abc" );
			TS_ASSERT( false );
		} catch ( utility::excn::Exception & e ) {
			//std::cout << "Error message:\n" << e.msg() << std::endl;
			std::string new_msg = e.msg();
			std::string target = "File: src/basic/resource_manager/ResourceManager.cc:";
			platform::Size start = new_msg.find(target) + target.size();
			for ( platform::Size ii = start; ii < new_msg.size(); ++ii ) {
				if ( new_msg[ ii ] >= '0' && new_msg[ ii ] <= '9' ) {
					new_msg[ ii ] = 'X';
				} else {
					break;
				}
			}

			std::string gold_standard_err_msg =
				"\n\nFile: src/basic/resource_manager/ResourceManager.cc:XXX\n"
				"Critical error: a resource may not depend on another resource that in turn depends on the first one. Infinite recursion is not allowed.\n"
				"Resources that were being created were:\n"
				"   1abc\n"
				"   2def\n"
				"   3ghi\n"
				"    ultimately requesting resource 1abc\n"
				"Cannot proceed from here.\n";

			TS_ASSERT_EQUALS( new_msg, gold_standard_err_msg );
		}

	}

};
