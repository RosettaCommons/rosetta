// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueenOverrideOutputter.cxxtest.hh
/// @brief  test suite for the StandardJobQueen testing when derived JQs intend to use
///         a limited subset of PoseOutputters.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/standard/StandardJobQueen.hh>

// Package headers
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/PoseInputSource.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

#include <protocols/moves/NullMover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

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

using namespace utility::tag;
using namespace protocols;
using namespace protocols::jd3;
using namespace protocols::jd3::pose_outputters;
using namespace protocols::jd3::standard;

//local options
namespace basic { namespace options { namespace OptionKeys {
basic::options::BooleanOptionKey const dummy_outputter_arg("dummy_outputter_arg");
basic::options::BooleanOptionKey const dummy_outputter("dummy_outputter"); // does the user want the dummy outputter?
}}}


// The derived job queen that overrides the standard set of outputters
class DerivedJobQueenOutOverride : public StandardJobQueen
{
public:
	DerivedJobQueenOutOverride() {}
	~DerivedJobQueenOutOverride() {}

	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP ,
		utility::options::OptionCollectionCOP ,
		utility::vector1< JobResultCOP > const &
	) override
	{
		core::pose::PoseOP pose( new core::pose::Pose );
		core::pose::make_pose_from_sequence( *pose, "ACDEFGH", core::chemical::FA_STANDARD );
		moves::MoverOP null_mover( new moves::NullMover );
		MoverAndPoseJobOP job( new MoverAndPoseJob );
		job->mover( null_mover );
		job->pose( pose );
		return job;
	}

	pose_outputters::PoseOutputterOP
	outputter( LarvalJobOP job ) {
		return pose_outputter_for_job( dynamic_cast< StandardInnerLarvalJob const & > ( *job->inner_job() ) );
	}
};

////
////

class DummyPoseOutputter : public PoseOutputter
{
public:
	std::map< std::string, core::pose::PoseOP > outputs;

	void determine_job_tag(
		utility::tag::TagCOP,
		utility::options::OptionCollection const &,
		InnerLarvalJob & job
	) const override {
		job.job_tag( job.input_tag() );
	}

	std::string
	outputter_for_job(
		utility::tag::TagCOP,
		utility::options::OptionCollection const &,
		InnerLarvalJob const &
	) const override {
		return "dummy";
	}

	bool job_has_already_completed( LarvalJob const &, utility::options::OptionCollection const & ) const override { return false; }

	void mark_job_as_having_started( LarvalJob const &, utility::options::OptionCollection const & ) const override {}

	void write_output_pose(
		LarvalJob const & job,
		utility::options::OptionCollection const &,
		utility::tag::TagCOP,
		core::pose::Pose const & pose
	) override {
		using namespace core::pose;
		outputs[ job.nstruct_suffixed_job_tag() ] = PoseOP( new Pose( pose ) );
	}

	void flush() override {}

	std::string
	class_key() const override { return "Dummy"; }

};

typedef utility::pointer::shared_ptr< DummyPoseOutputter > DummyPoseOutputterOP;

// An unregistered outputter
class DummyOutputterCreator : public PoseOutputterCreator
{
public:

	PoseOutputterOP create_outputter() const override { return PoseOutputterOP( new DummyPoseOutputter ); }
	std::string keyname() const override { return "Dummy"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override {
		AttributeList attrs;
		attrs + XMLSchemaAttribute( "dummy_attribute", xs_string, "A dummy attribute" );
		pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(), "testing 123", attrs );
	}
	void list_options_read( utility::options::OptionKeyList & read_options ) const override
	{
		read_options + basic::options::OptionKeys::dummy_outputter_arg;
	}

	bool outputter_specified_by_command_line() const override {
		return basic::options::option[ basic::options::OptionKeys::dummy_outputter ];
	}


};


// In the construtor for this class, the override functions for the PoseOutputter
// are called, so that the only allowed outputter for this class is the DummyOutputter.
class DerivedJobQueenOutOverride2 : public DerivedJobQueenOutOverride
{
public:
	DerivedJobQueenOutOverride2() {
		do_not_accept_all_pose_outputters_from_factory();
		allow_pose_outputter( PoseOutputterCreatorOP( new DummyOutputterCreator ));
	}
	~DerivedJobQueenOutOverride2() {}
};

// In the constructor for this class, an additional PoseOutputter is allowed, but
// all of the PoseOutputters from the PoseOutputterFactory are still allowed, also.
class DerivedJobQueenOutOverride3 : public DerivedJobQueenOutOverride
{
public:
	DerivedJobQueenOutOverride3() {
		allow_pose_outputter( PoseOutputterCreatorOP( new DummyOutputterCreator ));
	}
	~DerivedJobQueenOutOverride3() {}
};


class StandardJobQueenOverrideOutputterTests : public CxxTest::TestSuite
{
public:

	StandardJobQueenOverrideOutputterTests() : local_options_added_( false ) {}

	void setUp() {
		if ( ! local_options_added_ ) {
			using namespace basic::options;
			option.add( OptionKeys::dummy_outputter_arg, "" ).def(false);
			option.add( OptionKeys::dummy_outputter, "" ).def(false);
			local_options_added_ = true;
		}
	}

	void test_sjq_allow_only_non_factory_registered_outputter_read_cl()
	{
		// The DJQ2 should ignore the -s on the command line because the PDBPoseOutputter
		// should have been removed.
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueenOutOverride2 djq2;
		JobDigraphOP job_dag = djq2.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq2.determine_job_list( 1, 100 );
		LarvalJobOP larval_job1 = jobs.front();
		TS_ASSERT_EQUALS( jobs.size(), 1 );
		utility::vector1< protocols::jd3::JobResultCOP > results;
		JobOP job1 = djq2.mature_larval_job( larval_job1, results );
		CompletedJobOutput output1 = job1->run();
		djq2.completed_job_result( larval_job1, output1.second );
		PoseOutputterOP outputter = djq2.outputter( larval_job1 );
		DummyPoseOutputterOP dummy_outputter = utility::pointer::dynamic_pointer_cast< DummyPoseOutputter > ( outputter );
		TS_ASSERT( dummy_outputter );
		if ( ! dummy_outputter ) return;
		TS_ASSERT_EQUALS( dummy_outputter->outputs.count( "1ubq_0001" ), 1 );
	}

	void test_sjq_allow_only_non_factory_registered_outputter_read_jobdef_file()
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

		DerivedJobQueenOutOverride2 djq2;
		djq2.determine_preliminary_job_list_from_xml_file( jobdef_file );
		JobDigraphOP job_dag = djq2.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq2.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 5 );
		for ( auto larval_job : jobs ) {
			utility::vector1< protocols::jd3::JobResultCOP > results;
			JobOP job = djq2.mature_larval_job( larval_job, results );
			CompletedJobOutput output = job->run();
			djq2.completed_job_result( larval_job, output.second );
		}

		PoseOutputterOP outputter = djq2.outputter( jobs.front() );
		DummyPoseOutputterOP dummy_outputter = utility::pointer::dynamic_pointer_cast< DummyPoseOutputter > ( outputter );
		TS_ASSERT( dummy_outputter );
		if ( ! dummy_outputter ) return;
		for ( core::Size ii = 1; ii <= 5; ++ii ) {
			TS_ASSERT_EQUALS( dummy_outputter->outputs.count( "1ubq_000" + utility::to_string( ii ) ), 1 );
		}

	}


	void test_sjq_allow_only_non_factory_registered_outputter_read_jobdef_file_bad_outputter()
	{
		core_init_with_additional_options( "" );

		std::string jobdef_file =
			"<JobDefinitionFile>\n"
			" <Job nstruct=\"5\">\n"
			"  <Input>\n"
			"   <PDB filename=\"1ubq.pdb\"/>\n"
			"  </Input>\n"
			"  <Output>\n"
			"   <PDB filename=\"output.pdb\"/>\n"
			"  </Output>\n"
			" </Job>\n"
			"</JobDefinitionFile>\n";

		core_init(); // all options passed through job-definition file

		DerivedJobQueenOutOverride2 djq2;
		try {
			// This should fail because PDB shouldn't be an acceptible output type
			djq2.determine_preliminary_job_list_from_xml_file( jobdef_file );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string gold = "Job definition file failed to validate against the schema for this application\nUse the option -jd3::job_definition_schema <output filename> to output the schema to a file.\nError messages were: From line 7:\nError: Element \'PDB\': This element is not expected. Expected is ( Dummy ).\n\n 2:  <Job nstruct=\"5\">\n 3:   <Input>\n 4:    <PDB filename=\"1ubq.pdb\"/>\n 5:   </Input>\n 6:   <Output>\n 7:    <PDB filename=\"output.pdb\"/>\n 8:   </Output>\n 9:  </Job>\n10: </JobDefinitionFile>\n11: \n\nWarning messages were: \n";
			TS_ASSERT_EQUALS( e.msg(), gold );
		}
	}

	/////

	void test_sjq_add_unregistered_outputter_read_cl_get_default_outputter()
	{
		// The DJQ3 should prefer the PDBPoseOutputter since the user is not providing a request
		// to load the DummyOutputter.
		core_init_with_additional_options( "-s 1ubq.pdb" );
		DerivedJobQueenOutOverride3 djq3;
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );

		LarvalJobs jobs = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 1 );
		LarvalJobOP larval_job1 = jobs.front();

		PoseOutputterOP outputter = djq3.outputter( larval_job1 );
		PDBPoseOutputterOP pdb_outputter = utility::pointer::dynamic_pointer_cast< PDBPoseOutputter > ( outputter );
		TS_ASSERT( pdb_outputter );

	}

	void test_sjq_add_unregistered_outputter_read_cl_get_requested_dummy_outputter()
	{
		// The DJQ3 should prefer the PDBPoseOutputter since the user is not providing a request
		// to load the DummyOutputter.
		core_init_with_additional_options( "-s 1ubq.pdb -dummy_outputter" );
		DerivedJobQueenOutOverride3 djq3;
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );

		LarvalJobs jobs = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 1 );
		LarvalJobOP larval_job1 = jobs.front();

		PoseOutputterOP outputter = djq3.outputter( larval_job1 );
		DummyPoseOutputterOP dummy_outputter = utility::pointer::dynamic_pointer_cast< DummyPoseOutputter > ( outputter );
		TS_ASSERT( dummy_outputter );

	}

	void dont_test_sjq_add_unregistered_outputter_read_jobdef_file()
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

		DerivedJobQueenOutOverride3 djq3;
		djq3.determine_preliminary_job_list_from_xml_file( jobdef_file );
		JobDigraphOP job_dag = djq3.initial_job_dag();
		TS_ASSERT_EQUALS( job_dag->num_nodes(), 1 );
		LarvalJobs jobs = djq3.determine_job_list( 1, 100 );
		TS_ASSERT_EQUALS( jobs.size(), 5 );
		PoseInputSource const & input_source_1 = jobs.front()->inner_job()->input_source();
		TS_ASSERT_EQUALS( input_source_1.origin(), "Dummy" );
		TS_ASSERT_EQUALS( input_source_1.input_tag(), "dummy_from_tag" );
	}

	void dont_test_sjq_add_unregistered_outputter_read_jobdef_file2()
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

		DerivedJobQueenOutOverride3 djq3;
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
