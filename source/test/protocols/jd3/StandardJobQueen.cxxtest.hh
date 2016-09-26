// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueen.cxxtest.hh
/// @brief  test suite for the StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/StandardJobQueen.hh>

// Package headers
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

// basic headers
#include <basic/options/option.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <sstream>

using namespace utility::tag;
using namespace protocols::jd3;

//local options
namespace basic { namespace options { namespace OptionKeys {
basic::options::BooleanOptionKey const bool_arg1("bool_arg1");
basic::options::BooleanOptionKey const bool_arg2("bool_arg2");
basic::options::BooleanOptionKey const bool_arg3("bool_arg3");
basic::options::BooleanOptionKey const bool_arg_sjq_does_not_track("bool_arg_sjq_does_not_track");
basic::options::StringOptionKey const string_arg1("string_arg1");
basic::options::StringOptionKey const string_arg2("string_arg2");
basic::options::StringOptionKey const string_arg_w_default("string_arg_w_default");
basic::options::StringOptionKey const string_arg_sjq_does_not_track("string_arg_sjq_does_not_track");
namespace dummy_jq {
basic::options::BooleanOptionKey const bool_arg4("dummy_jq:bool_arg4");
}
utility::options::IntegerVectorOptionKey const intvect_arg1("intvect_arg1");
}}}//basic::options::OptionKeys

class DummyJobQueen : public ::StandardJobQueen
{
public:
	DummyJobQueen()
	{
		utility::options::OptionKeyList opts;

		add_options( opts );
		add_option( basic::options::OptionKeys::bool_arg1 );
		add_option( basic::options::OptionKeys::bool_arg2 );
		add_option( basic::options::OptionKeys::bool_arg3 );
		add_option( basic::options::OptionKeys::string_arg1 );
		add_option( basic::options::OptionKeys::string_arg2 );
		add_option( basic::options::OptionKeys::string_arg_w_default );
		add_option( basic::options::OptionKeys::dummy_jq::bool_arg4 );
		add_option( basic::options::OptionKeys::intvect_arg1 );
	}

	~DummyJobQueen() {}

	virtual
	void append_job_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen
	) const
	{
		if ( append_job_tag_subelements_ ) {
			append_job_tag_subelements_( job_definition_xsd, job_ct_gen );
		}
	}

	virtual
	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const
	{
		if ( append_common_tag_subelements_ ) {
			append_common_tag_subelements_( job_definition_xsd, ct_gen );
		}
	}


	virtual
	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< JobResultCOP > const &
	)
	{
		if ( complete_job_maturation_ ) {
			complete_job_maturation_( larval_job, job_options );
		}
		return protocols::jd3::JobOP();
	}

	//virtual bool has_job_completed( protocols::jd3::LarvalJobCOP job ) { return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job ); }
	virtual void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP /*job*/ ) {/*TEMP*/}

	virtual void note_job_completed( protocols::jd3::LarvalJobCOP /*job*/, protocols::jd3::JobStatus /*status*/ ) {}

	virtual void completed_job_result( protocols::jd3::LarvalJobCOP /*job*/, protocols::jd3::JobResultOP /*result*/ ) {}


public:

	// callbacks
	typedef boost::function< void ( protocols::jd3::LarvalJobCOP, utility::options::OptionCollectionCOP ) > CompleteJobMaturationCallback;
	CompleteJobMaturationCallback complete_job_maturation_;

	typedef boost::function< void ( utility::tag::XMLSchemaDefinition & xsd, utility::tag::XMLSchemaComplexTypeGenerator & ctgen ) > SubtagAppenderCallback;
	SubtagAppenderCallback append_job_tag_subelements_;
	SubtagAppenderCallback append_common_tag_subelements_;

};


class StandardJobQueenTests : public CxxTest::TestSuite
{
public:

	StandardJobQueenTests() : local_options_added_( false ) {}

	void setUp() {
		if ( ! local_options_added_ ) {
			using namespace basic::options;
			option.add( OptionKeys::bool_arg1, "" ).def(false);
			option.add( OptionKeys::bool_arg2, "" ).def(false);
			option.add( OptionKeys::bool_arg3, "" ).def(true);
			option.add( OptionKeys::bool_arg_sjq_does_not_track, "" ).def(false);
			option.add( OptionKeys::string_arg1, "" );
			option.add( OptionKeys::string_arg2, "" );
			option.add( OptionKeys::string_arg_sjq_does_not_track, "" );
			option.add( OptionKeys::string_arg_w_default, "" ).def("fiddlesticks");
			option.add( OptionKeys::dummy_jq::bool_arg4, "" ).def(false);
			option.add( OptionKeys::intvect_arg1, "" );

			local_options_added_ = true;
		}
	}

	void test_job_options_initialization() {
		core_init_with_additional_options( "-bool_arg1 -bool_arg_sjq_does_not_track -string_arg1 wakka_wakka_wakka -string_arg_sjq_does_not_track yippie -s 1ubq.pdb -intvect_arg1 1 2 3 4 5" );
		DummyJobQueen djq;
		djq.initial_job_dag(); // no need to hold the DAG returned by this func, but it must be called
		LarvalJobs jobs = djq.determine_job_list( 1, 1000 );
		TS_ASSERT( ! jobs.empty() );
		if ( jobs.empty() ) return;

		djq.complete_job_maturation_ = boost::bind( StandardJobQueenTests::callback_complete_larval_job_maturation1, _1, _2 );
		utility::vector1< JobResultOP > empty_vector;
		djq.mature_larval_job( jobs.front(), empty_vector ); // invokes callback_complete_larval_job_maturation1
	}

	static
	void
	callback_complete_larval_job_maturation1(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options_ptr
	)
	{
		using namespace basic::options::OptionKeys;
		TS_ASSERT_EQUALS( larval_job->inner_job()->n_input_sources(), 1 );
		TS_ASSERT_EQUALS( larval_job->inner_job()->input_source().origin(), pose_inputters::PDBPoseInputter::keyname() );
		TS_ASSERT_EQUALS( larval_job->inner_job()->input_source().input_tag(), "1ubq" );
		TS_ASSERT_EQUALS( larval_job->inner_job()->job_tag(), "1ubq" );

		TS_ASSERT( job_options_ptr );
		utility::options::OptionCollection const & job_options( * job_options_ptr );

		TS_ASSERT(   job_options[ bool_arg1 ].user() );
		TS_ASSERT( ! job_options[ bool_arg2 ].user() );
		TS_ASSERT( ! job_options[ bool_arg3 ].user() );

		TS_ASSERT(   job_options[ bool_arg1 ].active() );
		TS_ASSERT(   job_options[ bool_arg2 ].active() );
		TS_ASSERT(   job_options[ bool_arg3 ].active() );

		TS_ASSERT(   job_options[ bool_arg1 ] );
		TS_ASSERT( ! job_options[ bool_arg2 ] );
		TS_ASSERT(   job_options[ bool_arg3 ] );

		TS_ASSERT(   job_options[ string_arg1          ].user() );
		TS_ASSERT( ! job_options[ string_arg2          ].user() );
		TS_ASSERT( ! job_options[ string_arg_w_default ].user() );

		TS_ASSERT(   job_options[ string_arg1          ].active() );
		TS_ASSERT( ! job_options[ string_arg2          ].active() );
		TS_ASSERT(   job_options[ string_arg_w_default ].active() );

		TS_ASSERT(   job_options[ string_arg1 ]() == "wakka_wakka_wakka" );

		TS_ASSERT( ! job_options.has( bool_arg_sjq_does_not_track ) );
		TS_ASSERT( ! job_options.has( string_arg_sjq_does_not_track ) );

		TS_ASSERT(   job_options[ intvect_arg1         ].user() );
		TS_ASSERT(   job_options[ intvect_arg1         ].active() );
		utility::vector1< int > intvect_arg1_expected( 5 );
		intvect_arg1_expected[ 1 ] = 1;
		intvect_arg1_expected[ 2 ] = 2;
		intvect_arg1_expected[ 3 ] = 3;
		intvect_arg1_expected[ 4 ] = 4;
		intvect_arg1_expected[ 5 ] = 5;
	}

	bool tag_has_subtag_w_name( TagCOP tag, std::string const & tag_name, std::string const & name_attribute ) {
		for ( Tag::tags_t::const_iterator iter = tag->getTags().begin(); iter != tag->getTags().end(); ++iter ) {
			if ( (*iter)->getName() != tag_name ) continue;
			if ( (*iter)->hasOption( "name" ) ) {
				if ( (*iter)->getOption< std::string >( "name" ) == name_attribute ) {
					return true;
				}
			}
		}
		return false;
	}

	TagCOP
	subtag_w_name( TagCOP tag, std::string const & tag_name, std::string const & name_attribute ) {
		for ( Tag::tags_t::const_iterator iter = tag->getTags().begin(); iter != tag->getTags().end(); ++iter ) {
			if ( (*iter)->getName() != tag_name ) continue;
			if ( (*iter)->hasOption( "name" ) ) {
				if ( (*iter)->getOption< std::string >( "name" ) == name_attribute ) {
					return *iter;
				}
			}
		}
		return TagCOP();
	}

	void test_job_definition_file_xsd() {
		core_init();
		DummyJobQueen djq;
		djq.append_job_tag_subelements_ = boost::bind( StandardJobQueenTests::job_tag_xsd1, _1, _2 );
		std::string job_def_xsd = djq.job_definition_xsd();
		// now lets turn this into a Tag object and then make sure the Job tag has the definition it ought to

		//std::cout << "job def xsd:\n" << job_def_xsd << std::endl;

		TagCOP job_def_xsd_tag = Tag::create( job_def_xsd );
		TS_ASSERT_EQUALS( job_def_xsd_tag->getName(), std::string( "xs:schema" ) );

		TS_ASSERT( tag_has_subtag_w_name( job_def_xsd_tag, "xs:complexType", "job_def_Job_type" ));
		TS_ASSERT( tag_has_subtag_w_name( job_def_xsd_tag, "xs:complexType", "job_def_Options_type" ));
		TagCOP options_type_tag = subtag_w_name( job_def_xsd_tag, "xs:complexType", "job_def_Options_type" );
		if ( ! options_type_tag ) return;

		TS_ASSERT( options_type_tag->hasTag( "xs:all" ));
		TagCOP options_type_xs_all = options_type_tag->getTag( "xs:all" );
		if ( ! options_type_xs_all ) return;

		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "bool_arg1" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "bool_arg2" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "bool_arg3" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "string_arg1" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "string_arg2" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "string_arg_w_default" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "dummy_jq__bool_arg4" ));
		TS_ASSERT( tag_has_subtag_w_name( options_type_xs_all, "xs:element", "intvect_arg1" ));
	}

	static
	void job_tag_xsd1( utility::tag::XMLSchemaDefinition &, utility::tag::XMLSchemaComplexTypeGenerator & ctgen )
	{
		using namespace utility::tag;
		AttributeList attributes;
		attributes + XMLSchemaAttribute( "foo", xs_string ) + XMLSchemaAttribute( "bar", xs_integer );

		XMLSchemaSimpleSubelementList subelements;
		subelements.add_simple_subelement( "Hoo", attributes ).add_simple_subelement( "Ville", attributes );
		ctgen.add_ordered_subelement_set_as_pick_one( subelements );
	}

	void test_read_jobs_from_xml_file()
	{

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job>\n"
			"  <Input>\n"
			"   <PDB filename=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			"  <Options>\n"
			"   <bool_arg1/>\n"
			"   <string_arg1 value=\"wakka_wakka_wakka\"/>\n"
			"   <intvect_arg1 value=\"1 2 3 4 5\"/>\n"
			"  </Options>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DummyJobQueen djq;
		try {
			djq.determine_preliminary_job_list_from_xml_file( jobdef_file );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::cout << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		djq.initial_job_dag(); // no need to hold the DAG returned by this func, but it must be called

		LarvalJobs jobs = djq.determine_job_list( 1, 1000 );

		utility::vector1< JobResultOP > empty_vector;
		djq.complete_job_maturation_ = boost::bind( StandardJobQueenTests::callback_complete_larval_job_maturation1, _1, _2 );
		djq.mature_larval_job( jobs.front(), empty_vector ); // invokes callback_complete_larval_job_maturation1

	}


private:
	bool local_options_added_;

};
