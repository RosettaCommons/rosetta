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
#include <protocols/jd3/pose_outputters/ScoreFileOutputter.hh>

// Package headers
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


class ScoreFileOutputterTests : public CxxTest::TestSuite
{
public:

	void setUp() {
	}

	void test_output_path_handling_1() {
		core_init_with_additional_options( "" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "dummy" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "score.sc" );

	}

	void test_output_path_handling_2() {
		core_init_with_additional_options( "-out:file:scorefile dummy.sc" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "dummy.sc" );

	}

	void test_output_path_handling_3() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "dummy.scores" );

	}

	void test_output_path_handling_4a() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:prefix 2def_" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "dummy.scores" );

	}

	void test_output_path_handling_4b() {
		core_init_with_additional_options( "-out:prefix 2def_" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "2def_score.sc" );

	}


	void test_output_path_handling_5a() {
		core_init_with_additional_options( "-out:suffix _varcst" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "score_varcst.sc" );

	}

	void test_output_path_handling_5b() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:suffix _varcst" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "dummy.scores" );

	}

	void test_output_path_handling_6() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path some_dir1" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "some_dir1/dummy.scores" );

	}

	void test_output_path_handling_7() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path:all some_dir2" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "some_dir2/dummy.scores" );

	}

	void test_output_path_handling_8() {
		// no name given for out::file::scorefile
		core_init_with_additional_options( "-out:path:all some_dir2" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "some_dir2/score.sc" );

	}

	void test_output_path_handling_priority_1() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path:score some_dir3 -out:path:all some_dir4" );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );
		utility::tag::TagCOP null_tag;

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( null_tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "some_dir3/dummy.scores" );

	}

	void test_output_path_handling_priority_2() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path:score some_dir3 -out:path:all some_dir4" );

		std::string jobdef_snippet = "   <ScoreFile filename=\"score1.sc\"/>\n";
		utility::tag::TagCOP tag = utility::tag::Tag::create( jobdef_snippet );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "some_dir3/score1.sc" );

	}

	void test_output_path_handling_priority_3() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path:score some_dir3 -out:path:all some_dir4" );

		std::string jobdef_snippet = "   <ScoreFile filename=\"./score1.sc\"/>\n";
		utility::tag::TagCOP tag = utility::tag::Tag::create( jobdef_snippet );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "./score1.sc" );

	}

	void test_output_path_handling_priority_4() {
		core_init_with_additional_options( "-out:file:scorefile dummy.scores -out:path:score some_dir3 -out:path:all some_dir4" );

		std::string jobdef_snippet = "   <ScoreFile path=\"tagdir\"/>\n";
		utility::tag::TagCOP tag = utility::tag::Tag::create( jobdef_snippet );

		ScoreFileOutputter outputter;
		utility::options::OptionKeyList pdb_outputter_options;
		ScoreFileOutputter::list_options_read( pdb_outputter_options );

		auto dummy_input_source = utility::pointer::make_shared< PoseInputSource >( "PDB" );
		dummy_input_source->input_tag( "1abc" );
		auto inner_job = utility::pointer::make_shared< InnerLarvalJob >( 1 );
		inner_job->input_source( dummy_input_source );

		utility::options::OptionCollectionOP job_options = basic::options::read_subset_of_global_option_collection( pdb_outputter_options );

		JobOutputIndex output_index;
		// initialize the output index

		utility::file::FileName output_name = outputter.filename_for_job( tag, *job_options, *inner_job );
		TS_ASSERT_EQUALS( output_name(), "tagdir/dummy.scores" );

	}

};
