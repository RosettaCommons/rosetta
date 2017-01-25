// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueenOverrideInputter.cxxtest.hh
/// @brief  test suite for the StandardJobQueen testing when derived JQs intend to use
///         a limited subset of PoseInputters.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/standard/StandardJobQueen.hh>

// Package headers
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>
#include <protocols/jd3/pose_inputters/pose_inputter_schemas.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>
#include <protocols/jd3/pose_inputters/PoseInputterFactory.hh>
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <sstream>
#include <algorithm>

using namespace utility::tag;
using namespace protocols::jd3;
using namespace protocols::jd3::pose_inputters;
using namespace protocols::jd3::standard;

//local options
namespace basic { namespace options { namespace OptionKeys {
basic::options::BooleanOptionKey const dummy_inputter_arg("dummy_inputter_arg");
}}}


// The derived job queen that overrides the standard set of inputters
class DerivedJobQueen : public StandardJobQueen
{
public:
	DerivedJobQueen() {}
	~DerivedJobQueen() {}

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP ,
		utility::options::OptionCollectionCOP ,
		utility::vector1< JobResultCOP > const &
	) override
	{
		return protocols::jd3::JobOP();
	}
};

////
////

class DummyPoseInputter : public PoseInputter
{
	bool job_available_on_command_line() const override { return true; }
	PoseInputSources pose_input_sources_from_command_line() override {
		PoseInputSources one_source;
		PoseInputSourceOP pis( new PoseInputSource( "Dummy" ) );
		pis->input_tag( "dummy_from_cl" );
		one_source.push_back( pis );
		return one_source;
	}

	PoseInputSources
	pose_input_sources_from_tag(
		utility::options::OptionCollection const &,
		utility::tag::TagCOP
	) override {
		PoseInputSources one_source;
		PoseInputSourceOP pis( new PoseInputSource( "Dummy" ) );
		pis->input_tag( "dummy_from_tag" );
		one_source.push_back( pis );
		return one_source;
	}

	core::pose::PoseOP
	pose_from_input_source(
		PoseInputSource const &,
		utility::options::OptionCollection const &,
		utility::tag::TagCOP
	) override {
		return core::pose::PoseOP();
	}
};

// An unregistered inputter
class DummyInputterCreator : public PoseInputterCreator
{
public:
	PoseInputterOP create_inputter() const { return PoseInputterOP( new DummyPoseInputter ); }
	std::string keyname() const { return "Dummy"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
		AttributeList attrs;
		attrs + XMLSchemaAttribute( "dummy_attribute", xs_string, "A dummy attribute" );
		pose_inputter_xsd_type_definition_w_attributes( xsd, keyname(), "testing 123", attrs );
	}
	void list_options_read( utility::options::OptionKeyList & read_options ) const
	{
		read_options + basic::options::OptionKeys::dummy_inputter_arg;
	}
};


// In the construtor for this class, the override functions for the PoseInputter
// are called, so that the only allowed inputter for this class is the DummyInputter.
class DerivedJobQueen2 : public DerivedJobQueen
{
public:
	DerivedJobQueen2() {
		do_not_accept_all_pose_inputters_from_factory();
		allow_pose_inputter( PoseInputterCreatorOP( new DummyInputterCreator ));
	}
	~DerivedJobQueen2() {}
};

// In the constructor for this class, an additional PoseInputter is allowed, but
// all of the PoseInputters from the PoseInputterFactory are still allowed, also.
class DerivedJobQueen3 : public DerivedJobQueen
{
public:
	DerivedJobQueen3() {
		allow_pose_inputter( PoseInputterCreatorOP( new DummyInputterCreator ));
	}
	~DerivedJobQueen3() {}
};


class StandardJobQueenOverrideInputterTests : public CxxTest::TestSuite
{
public:

	StandardJobQueenOverrideInputterTests() : local_options_added_( false ) {}

	void setUp() {
		if ( ! local_options_added_ ) {
			using namespace basic::options;
			option.add( OptionKeys::dummy_inputter_arg, "" ).def(false);
			local_options_added_ = true;
		}
	}


	template < class It, class Pr, class T >
	T
	find_if(
		It begin,
		It end,
		Pr predicate,
		T not_found_return_val
	)
	{
		It it = std::find_if( begin, end, predicate );
		if ( it == end ) {
			return not_found_return_val;
		} else {
			return *it;
		}
	}

	template < class Container, class Pr, class T >
	T
	find_if(
		Container const & container,
		Pr predicate,
		T not_found_return_val
	)
	{
		return find_if( container.begin(), container.end(), predicate, not_found_return_val );
	}

	TagCOP
	find_inputter_tag(
		TagCOP tag,
		std::string const & inputter_name
	){
		std::string ct_name = PoseInputterFactory::complex_type_name_for_pose_inputter( inputter_name );
		return find_if(
			tag->getTags( "xs:complexType" ),
			[ &ct_name ] ( TagCOP const & subtag ) { return subtag->getOption< std::string > ( "name", "" ) == ct_name; },
			TagCOP() );
	}

	TagCOP
	find_option_tag(
		TagCOP tag, // should be the top level of the XSD
		utility::options::OptionKey const & option_key
	) {
		std::string decolonized_option_name = basic::options::replace_option_namespace_colons_with_underscores( option_key );

		std::string opt_ct_name = StandardJobQueen::job_def_complex_type_name( "Options" );
		TagCOP option_ct = find_if(
			tag->getTags( "xs:complexType" ),
			[ &opt_ct_name ] ( TagCOP const & subtag ) { return subtag->getOption< std::string > ( "name", "" ) == opt_ct_name; },
			TagCOP() );
		if ( ! option_ct ) { return option_ct; }

		return find_if(
			option_ct->getTag( "xs:all" )->getTags( "xs:element" ),
			[ &decolonized_option_name ] ( TagCOP const & subtag ) { return subtag->getOption< std::string > ( "name", "" ) == decolonized_option_name; },
			TagCOP() );
	}

	void test_sjq_allow_only_non_factory_registered_inputter_read_cl()
	{
		// The DJQ2 should ignore the -s on the command line because the PDBPoseInputter
		// should have been removed.
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueen2 djq2;
		JobDigraphOP job_dag = djq2.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq2.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 1 );
		PoseInputSource const & input_source_1 = jobs.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "dummy_from_cl" );
	}

	void test_sjq_allow_only_non_factory_registered_inputter_xsd_correct()
	{
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueen2 djq2;
		std::string xsd = djq2.job_definition_xsd();
		TagOP xsd_tag = Tag::create( xsd );
		TagCOP dummy_inputter_ct_tag = find_inputter_tag( xsd_tag, "Dummy" );
		TS_ASSERT( dummy_inputter_ct_tag );
		TagCOP pdb_inputter_ct_tag = find_inputter_tag( xsd_tag, "PDB" );
		TS_ASSERT( ! pdb_inputter_ct_tag );

		using namespace basic::options::OptionKeys;

		TagCOP dummy_opt_tag = find_option_tag( xsd_tag, dummy_inputter_arg );
		TS_ASSERT( dummy_opt_tag );

		TagCOP in_file_nco_opt_tag = find_option_tag( xsd_tag, in::file::new_chain_order );
		TS_ASSERT( ! in_file_nco_opt_tag ); // As long as the PDB inputter is not used,
		// this option key should not be found in the list of acceptible options.

	}

	void test_sjq_allow_only_non_factory_registered_inputter_read_jobdef_file()
	{
		core_init_with_additional_options( "" );

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <Dummy dummy_attribute=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DerivedJobQueen2 djq2;
		djq2.determine_preliminary_job_list_from_xml_file( jobdef_file );
		JobDigraphOP job_dag = djq2.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq2.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 5 );
		PoseInputSource const & input_source_1 = jobs.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "dummy_from_tag" );

	}

	void test_sjq_allow_only_non_factory_registered_inputter_read_jobdef_file_bad_inputter()
	{
		core_init_with_additional_options( "" );

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <PDB filename=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DerivedJobQueen2 djq2;
		try {
			// This should fail because PDB shouldn't be an acceptible input type
			djq2.determine_preliminary_job_list_from_xml_file( jobdef_file );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string gold = "Job definition file failed to validate against the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output the schema to a file.\nError messages were: From line 4:\nError: Element \'PDB\': This element is not expected. Expected is ( Dummy ).\n\n1: <JobDefinitionFile>\n2:  <Job nstruct=\"5\">\n3:   <Input>\n4:    <PDB filename=\"1ubq.pdb\"/>\n5:   </Input>\n6:  </Job>\n7: </JobDefinitionFile>\n8: \n\nWarning messages were: \n";
			TS_ASSERT_EQUALS( e.msg(), gold );
		}
	}

	/////

	void test_sjq_add_unregistered_inputter_read_cl()
	{
		// The DJQ3 should not ignore the -s on the command line
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueen3 djq3;
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 2 );

		LarvalJobs jobs1 = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs1.size(), 1 );
		PoseInputSource const & input_source_1 = jobs1.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "PDB" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "1ubq" );

		LarvalJobs jobs2 = djq3.determine_job_list( 2, 100 );
		TS_ASSERT_EQUALS( jobs2.size(), 1 );
		PoseInputSource const & input_source_2 = jobs2.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_2.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_2.input_tag(), "dummy_from_cl" );
	}

	void test_sjq_add_unregistered_inputter_xsd_correct()
	{
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueen3 djq3;
		std::string xsd = djq3.job_definition_xsd();
		//std::cout << xsd << "\n--------\n" << std::endl;
		TagOP xsd_tag = Tag::create( xsd );
		TagCOP dummy_inputter_ct_tag = find_inputter_tag( xsd_tag, "Dummy" );
		TS_ASSERT( dummy_inputter_ct_tag );
		TagCOP pdb_inputter_ct_tag = find_inputter_tag( xsd_tag, "PDB" );
		TS_ASSERT( pdb_inputter_ct_tag );

		using namespace basic::options::OptionKeys;

		TagCOP dummy_opt_tag = find_option_tag( xsd_tag, dummy_inputter_arg );
		TS_ASSERT( dummy_opt_tag );

		TagCOP in_file_nco_opt_tag = find_option_tag( xsd_tag, in::file::new_chain_order );
		TS_ASSERT( in_file_nco_opt_tag );
	}


	void test_sjq_add_unregistered_inputter_read_jobdef_file()
	{
		core_init_with_additional_options( "" );

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <Dummy dummy_attribute=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DerivedJobQueen3 djq3;
		djq3.determine_preliminary_job_list_from_xml_file( jobdef_file );
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 5 );
		PoseInputSource const & input_source_1 = jobs.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "dummy_from_tag" );
	}

	void test_sjq_add_unregistered_inputter_read_jobdef_file2()
	{
		core_init_with_additional_options( "" );

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <Dummy dummy_attribute=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			" </Job>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <PDB filename=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DerivedJobQueen3 djq3;
		djq3.determine_preliminary_job_list_from_xml_file( jobdef_file );
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 2 );
		LarvalJobs jobs1 = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs1.size(), 5 );
		PoseInputSource const & input_source_1 = jobs1.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "dummy_from_tag" );

		LarvalJobs jobs2 = djq3.determine_job_list( 2, 100 );
		TS_ASSERT_EQUALS( jobs2.size(), 5 );
		PoseInputSource const & input_source_2 = jobs2.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_2.origin(), "PDB" );
		TS_ASSERT_EQUALS( input_source_2.input_tag(), "1ubq" );

	}


private:
	bool local_options_added_;

};
