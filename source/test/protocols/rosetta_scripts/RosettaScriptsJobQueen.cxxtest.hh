// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/rosetta_scripts/RosettaScriptsJobQueen.cxxtest.hh
/// @brief  test suite for the RosettaScriptsJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/rosetta_scripts/RosettaScriptsJobQueen.hh>

// Package headers
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/LarvalJob.hh>
//#include <protocols/jd3/InnerLarvalJob.hh>
//#include <protocols/jd3/pose_inputters/PoseInputSource.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>
#include <protocols/jd3/deallocation/ResourceDeallocationMessage.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/moves/mover_schemas.hh>


// basic headers
#include <basic/options/option.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <basic/resource_manager/ResourceLoaderRegistrator.hh>
#include <basic/resource_manager/loader_schemas.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Boost headers
//#include <boost/bind.hpp>
//#include <boost/function.hpp>

// C++ headers
#include <sstream>

using namespace utility::tag;
using namespace protocols::jd3;
using namespace protocols::jd3::deallocation;
using namespace protocols::jd3::standard;
using namespace protocols::jd3::pose_inputters;
using namespace protocols::rosetta_scripts;

class DummyResource : public utility::pointer::ReferenceCount
{
public:
	DummyResource( int value ) : value_( value ) {}
	int value_;
};

class DummyResourceLoader : public basic::resource_manager::ResourceLoader
{
public:
	DummyResourceLoader() {}

	static
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
		utility::tag::AttributeList attlist;
		basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd, "Dummy", "(desc)", attlist );
	}

	basic::resource_manager::ResourceCOP
	create_resource(
		basic::resource_manager::ResourceManager &,
		utility::tag::TagCOP,
		std::string const &,
		std::istream &
	) const override
	{
		return std::make_shared< DummyResource >( ++count_resources );
	}

	static int count_resources;

};

int
DummyResourceLoader::count_resources( 0 );


class DummyResourceLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
{
public:
	/// @brief Return a up-casted owning pointer (ResourceLoaderOP) to the resource loader.
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const override {
		return std::make_shared< DummyResourceLoader >();
	}

	/// @brief Return the string identifier for the associated ResourceLoader (LoopsFile).
	std::string loader_type() const override {
		return "Dummy";
	}

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override
	{
		DummyResourceLoader::provide_xml_schema( xsd );
	}

};

class NullMover1ABC : public protocols::moves::NullMover
{
public:
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) override {
		//std::cout << "parse_my_tag for NullMover1ABC" << std::endl;
		TS_ASSERT( data.has_resource( "1abc" ) );
	}
};

class NullMover2DEF : public protocols::moves::NullMover
{
public:
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap & data,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) override {
		TS_ASSERT( data.has_resource( "2def" ) );
	}
};

class NullMover1ABCCreator : public protocols::moves::MoverCreator
{
	protocols::moves::MoverOP create_mover() const override { return std::make_shared< NullMover1ABC >(); }
	std::string keyname() const override { return "Null1ABC"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override {
		utility::tag::AttributeList attlist;
		protocols::moves::xsd_type_definition_w_attributes( xsd, "Null1ABC", "dummy", attlist );
	}
};


class NullMover2DEFCreator : public protocols::moves::MoverCreator
{
	protocols::moves::MoverOP create_mover() const override { return std::make_shared< NullMover2DEF >(); }
	std::string keyname() const override { return "Null2DEF"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override {
		utility::tag::AttributeList attlist;
		protocols::moves::xsd_type_definition_w_attributes( xsd, "Null2DEF", "dummy", attlist );
	}
};


class RosettaScriptsJobQueenTests : public CxxTest::TestSuite
{
private:
	bool registered_null_movers_ = { false };
public:

	void setUp() {
		if ( ! registered_null_movers_ ) {
			registered_null_movers_ = true;
			protocols::moves::MoverRegistrator< NullMover1ABCCreator > r1;
			protocols::moves::MoverRegistrator< NullMover2DEFCreator > r2;
			basic::resource_manager::ResourceLoaderRegistrator< DummyResourceLoaderCreator > r3;
		}
		core_init();
	}


	void test_read_resources_from_xml_file()
	{

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"2\">\n"
			"  <Input>\n"
			"   <PDB filename=\"protocols/rosetta_scripts/1abc.pdb\"/>\n"
			"  </Input>\n"
			"  <Options>\n"
			"   <parser__protocol value=\"protocols/rosetta_scripts/example_resource_using_script.xml\"/>\n"
			"   <parser__script_vars value=\"mover=Null1ABC native=1abc\"/>\n"
			"  </Options>\n"
			" </Job>\n"
			" <Job nstruct=\"2\">\n"
			"  <Input>\n"
			"   <PDB filename=\"protocols/rosetta_scripts/2def.pdb\"/>\n"
			"  </Input>\n"
			"  <Options>\n"
			"   <parser__protocol value=\"protocols/rosetta_scripts/example_resource_using_script.xml\"/>\n"
			"   <parser__script_vars value=\"mover=Null2DEF native=2def\"/>\n"
			"  </Options>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

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

		RosettaScriptsJobQueen jq1;
		jq1.resource_manager()->read_resources_from_xml( "(unit test)", resource_def_file );
		TS_ASSERT_EQUALS( jq1.resource_manager()->n_resources_declared(), 2 );
		try {
			jq1.determine_preliminary_job_list_from_xml_file( jobdef_file );
		} catch (utility::excn::Exception & e ) {
			std::cout << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		JobDigraphOP dag = jq1.initial_job_dag(); // no need to hold the DAG returned by this func, but it must be called
		TS_ASSERT_EQUALS( dag->num_nodes(), 2 );

		LarvalJobs jobs1 = jq1.determine_job_list( 1, 1000 );
		TS_ASSERT_EQUALS( jobs1.size(), 2 );

		LarvalJobs jobs2 = jq1.determine_job_list( 2, 1000 );
		TS_ASSERT_EQUALS( jobs2.size(), 2 );

		utility::vector1< JobResultOP > empty_vector;
		for ( auto const & larval_job : jobs1 ) {
			jq1.mature_larval_job( larval_job, empty_vector );
		}
		TS_ASSERT_EQUALS( jq1.resource_manager()->n_resources_in_memory(), 1 );

		for ( auto const & larval_job : jobs1 ) {
			jq1.note_job_completed( larval_job, jd3_job_status_success, 1 );
		}

		// OK - now, if we ask the JQ1 for her deallocation messages, she should say
		// that the 1abc resource needs to be deallocated.
		bool deallocated_1abc_pose( false ), deallocated_1abc_resource( false );
		std::list< DeallocationMessageOP > msgs = jq1.deallocation_messages();
		DeallocationMessageOP deallocate_1abc_resource_message;
		TS_ASSERT_EQUALS( msgs.size(), 2 );
		for ( auto const & msg : msgs ) {
			if ( msg->deallocation_type() == input_pose_deallocation_msg ) {
				deallocated_1abc_pose = true;
			} else if ( msg->deallocation_type() == resource_deallocation_msg ) {
				auto rdm = utility::pointer::dynamic_pointer_cast< ResourceDeallocationMessage >( msg );
				TS_ASSERT_EQUALS( rdm->resource_name(), "1abc" );
				deallocated_1abc_resource = true;
				deallocate_1abc_resource_message = msg;
			} else {
				// There should be
				TS_ASSERT( false );
			}
		}
		TS_ASSERT( deallocated_1abc_pose );
		TS_ASSERT( deallocated_1abc_resource );
		TS_ASSERT( deallocate_1abc_resource_message );
		TS_ASSERT_EQUALS( jq1.resource_manager()->n_resources_in_memory(), 0 );

		for ( auto const & larval_job : jobs2 ) {
			jq1.mature_larval_job( larval_job, empty_vector );
		}
		TS_ASSERT_EQUALS( jq1.resource_manager()->n_resources_in_memory(), 1 );

		for ( auto const & larval_job : jobs2 ) {
			jq1.note_job_completed( larval_job, jd3_job_status_success, 1 );
		}

		// OK - now, if we ask the JQ for her deallocation messages, she should say
		// that the 2def resource needs to be deallocated.
		// Currently, the logic in the StandardJobQueen is to not deallocate the
		// input pose for the last one (it's difficult for the base class to know
		// that no subsequent jobs will be asked for by the derived job queen until
		// the derived job queen goes on to the next input pose, and that cannot
		// happen for the very last input pose) so we don't here expect to see
		// that the 2def pose will be deleted.
		bool deallocated_2def_pose( false ), deallocated_2def_resource( false );
		std::list< DeallocationMessageOP > msgs2 = jq1.deallocation_messages();
		DeallocationMessageOP deallocate_2def_resource_message;
		TS_ASSERT_EQUALS( msgs2.size(), 1 );
		for ( auto const & msg : msgs2 ) {
			if ( msg->deallocation_type() == input_pose_deallocation_msg ) {
				//deallocated_1abc_pose = true;
			} else if ( msg->deallocation_type() == resource_deallocation_msg ) {
				auto rdm = utility::pointer::dynamic_pointer_cast< ResourceDeallocationMessage >( msg );
				TS_ASSERT_EQUALS( rdm->resource_name(), "2def" );
				deallocated_2def_resource = true;
				deallocate_2def_resource_message = msg;
			} else {
				// There should be
				TS_ASSERT( false );
			}
		}
		TS_ASSERT( ! deallocated_2def_pose );
		TS_ASSERT( deallocated_2def_resource );
		TS_ASSERT_EQUALS( jq1.resource_manager()->n_resources_in_memory(), 0 );

		///////////////////////////////////////////////////////////////////////////////////////////
		// OK! Let's create a second job queen that we will pretend is on a remote node
		// so we'll ask this job queen to mature a larval job for each of the 1abc and
		// 2def preliminary job nodes.
		//
		RosettaScriptsJobQueen jq2;
		jq2.resource_manager()->read_resources_from_xml( "(unit test)", resource_def_file );
		try {
			jq2.determine_preliminary_job_list_from_xml_file( jobdef_file );
		} catch (utility::excn::Exception & e ) {
			std::cout << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT_EQUALS( jq2.resource_manager()->n_resources_declared(), 2 );
		jq2.mature_larval_job( jobs1.front(), empty_vector );
		TS_ASSERT_EQUALS( jq2.resource_manager()->n_resources_in_memory(), 1 );

		jq2.mature_larval_job( jobs2.front(), empty_vector );
		TS_ASSERT_EQUALS( jq2.resource_manager()->n_resources_in_memory(), 2 );

		jq2.process_deallocation_message( deallocate_1abc_resource_message );
		TS_ASSERT_EQUALS( jq2.resource_manager()->n_resources_in_memory(), 1 );

		jq2.process_deallocation_message( deallocate_2def_resource_message );
		TS_ASSERT_EQUALS( jq2.resource_manager()->n_resources_in_memory(), 0 );


	}

};
