// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/pose_outputters/SilentFilePoseOutputter.cxxtest.hh
/// @brief  test suite for the class responsible for writing Poses to disk in silent file format
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputter.hh>

// Package headers
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.hh>

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

// C++ headers
#include <sstream>

using namespace utility::tag;
using namespace protocols::jd3;
using namespace protocols::jd3::pose_inputters;
using namespace protocols::jd3::pose_outputters;


class SilentFilePoseOutputterTests : public CxxTest::TestSuite
{
public:

	SilentFilePoseOutputterTests() {}

	void setUp() {
	}

	void test_output_path_handling_1a() {
		core_init_with_additional_options( "-out::file::silent dummy.out" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));
		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "dummy.out" );

	}

	void test_output_path_handling_1b() {
		// look for default ".out" extension
		core_init_with_additional_options( "-out::file::silent dummy" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));
		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "dummy.out" );

		//PoseOutputSpecificationOP output_spec = outputter.output_pdb_name( job, output_index, *job_options, null_tag );
		//TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "dummy.out" );

	}

	void test_output_path_handling_1c() {
		// if another exension is provided, it should be respected
		core_init_with_additional_options( "-out::file::silent dummy.silent" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "dummy.silent" );

	}

	void test_output_path_handling_2() {
		// output path for PDBs is going to be ignored
		core_init_with_additional_options( "-out::file::silent dummy -out::path::pdb some_dir" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;


		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "dummy.out" );

	}

	void test_output_path_handling_3() {
		core_init_with_additional_options( "-out::file::silent dummy -out::path::all some_dir2" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "some_dir2/dummy.out" );

	}

	void test_output_path_handling_4() {
		core_init_with_additional_options( "-out::file::silent dummy -out::path some_dir3" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "some_dir3/dummy.out" );

	}

	void test_output_path_handling_5() {
		core_init_with_additional_options( "-out::file::silent some_dir5/dummy" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "some_dir5/dummy.out" );

	}

	void test_output_path_handling_6() {
		core_init_with_additional_options( "-out::file::silent some_dir6/dummy -out:path root_destination" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "root_destination/some_dir6/dummy.out" );

	}

	void test_output_path_handling_7() {
		core_init_with_additional_options( "-out::file::silent /some_absolute_path/dummy -out:path root_destination" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "/some_absolute_path/dummy.out" );

	}


	void test_output_path_handling_priority_1() {
		core_init_with_additional_options( "-out::file::silent dummy -out::path some_dir3 -out::path::all some_dir4" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, null_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "some_dir4/dummy.out" );

	}

	void test_output_path_handling_with_tag() {
		core_init_with_additional_options( "-out::file::silent dummy -out::path some_dir3" );

		auto outputter_tag = utility::tag::Tag::create( "<SilentFile filename=\"foo.out\"/>" );

		SilentFilePoseOutputter outputter;
		utility::options::OptionKeyList sf_outputter_options;
		SilentFilePoseOutputter::list_options_read( sf_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );
		LarvalJob job( inner_job, 1, 1 );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( sf_outputter_options );

		JobOutputIndex output_index;

		PoseOutputSpecificationOP output_spec = outputter.create_output_specification( job, output_index, *job_options, outputter_tag );
		auto silent_output_spec( utility::pointer::dynamic_pointer_cast< SilentFilePoseOutputSpecification >( output_spec ));

		TS_ASSERT_EQUALS( silent_output_spec->out_fname(), "some_dir3/foo.out" );

	}

};
