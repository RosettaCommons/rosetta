// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/pose_outputters/PDBPoseOutputter.cxxtest.hh
/// @brief  test suite for the class responsible for writing Poses to disk in PDB format.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
#include <protocols/jd3/pose_outputters/mmTFPoseOutputter.hh>

// Package headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.hh>
#include <protocols/jd3/job_results/PoseJobResult.hh>
#include <protocols/jd3/pose_outputters/mmTFPoseOutputSpecification.hh>
#include <core/import_pose/import_pose.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <utility/file/file_sys_util.hh>

// C++ headers
#include <sstream>

using namespace utility::tag;
using namespace protocols::jd3;
using namespace protocols::jd3::pose_inputters;
using namespace protocols::jd3::pose_outputters;
using namespace protocols::jd3::job_results;


class PDBPoseOutputterTests : public CxxTest::TestSuite
{
public:

	PDBPoseOutputterTests() {}

	void setUp() {
	}

	void test_output_path_handling_1() {
		core_init_with_additional_options( "" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1, 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "dummy_0001.pdb" );

	}

	void test_output_path_handling_2() {
		core_init_with_additional_options( "-out::path::pdb some_dir" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1,1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir/dummy_0001.pdb" );

	}

	void test_output_path_handling_3() {
		core_init_with_additional_options( "-out::path::all some_dir2" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1,1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir2/dummy_0001.pdb" );

	}

	void test_output_path_handling_4() {
		core_init_with_additional_options( "-out::path some_dir3" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1, 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir3/dummy_0001.pdb" );

	}

	void test_output_path_handling_prefix_suffix() {
		core_init_with_additional_options( "-out::path some_dir3 -out:prefix pre_ -out:suffix _end" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1, 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		outputter.determine_job_tag(null_tag, *job_options, *inner_job);
		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir3/pre_dummy_end_0001.pdb" );

	}

	void test_output_path_tag_prefix_suffix(){
		core_init_with_additional_options( "-out::path some_dir3" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1, 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		std::string tag_str =
			"<Output>"
			"\t<PDB prefix=\"pre_\" suffix=\"_end\"/>"
			"</Output>";

		utility::tag::TagOP out_tag = utility::tag::Tag::create(tag_str);
		utility::tag::TagCOP null_tag;
		JobOutputIndex output_index;
		// initialize the output index

		outputter.determine_job_tag(out_tag, *job_options, *inner_job);

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir3/pre_dummy_end_0001.pdb" );

	}

	void test_output_path_tag_dir(){
		core_init_with_additional_options( "" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1, 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );


		std::string tag_str =
			"<Output>"
			"\t<PDB prefix=\"pre_\" suffix=\"_end\" filename_pattern=\"some_dir3/additional_$\"/>"
			"</Output>";

		utility::tag::TagOP out_tag = utility::tag::Tag::create(tag_str);
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		outputter.determine_job_tag(out_tag, *job_options, *inner_job);
		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir3/pre_additional_dummy_end_0001.pdb" );
	}

	void test_output_path_handling_priority_1() {
		core_init_with_additional_options( "-out::path some_dir3 -out::path::pdb some_dir4" );

		PDBPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		PDBPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1,1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir4/dummy_0001.pdb" );

	}
	void test_mmTF_output() {
		core_init_with_additional_options( "-out::path some_dir3" );

		mmTFPoseOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		mmTFPoseOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1,1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		TS_ASSERT_EQUALS( outputter.class_key(), "mmTF");

		std::string output_name = outputter.output_name( job, output_index, *job_options, null_tag );
		TS_ASSERT_EQUALS( output_name, "some_dir3/dummy_0001.mmtf" );

		//Test Actual Output
		core::pose::PoseOP pdb_pose = core::import_pose::pose_from_file( "core/io/1QYS.pdb", false , core::import_pose::PDB_file);
		PoseJobResult result = PoseJobResult(pdb_pose);
		mmTFPoseOutputSpecification spec = mmTFPoseOutputSpecification();
		spec.out_fname(output_name);

		TS_ASSERT_THROWS_NOTHING(outputter.write_output(spec, result));
		TS_ASSERT( utility::file::file_exists(output_name));


	}

};
