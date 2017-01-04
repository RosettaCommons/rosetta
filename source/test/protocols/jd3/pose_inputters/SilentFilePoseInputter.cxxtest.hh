// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/SilentFilePoseInputter.cxxtest.hh
/// @brief  test suite for protocols::jd3::SilentFilePoseInputer
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/jd3/pose_inputters/SilentFilePoseInputter.hh>
#include <protocols/jd3/PoseInputSource.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>

template <typename T>
typename T::mapped_type get(T const& map, typename T::key_type const& key)
{
	typename T::const_iterator iter(map.find(key));
	return iter != map.end() ? iter->second : typename T::mapped_type();
}

using protocols::jd3::PoseInputSources;
using protocols::jd3::PoseInputSource;
using protocols::jd3::PoseInputSourceOP;
using protocols::jd3::pose_inputters::SilentFilePoseInputter;

class SilentFilePoseInputterTests : public CxxTest::TestSuite {
public:

	void setUp() {
	}

	void test_read_in_file_silent_flag() {
		protocols_init_with_additional_options( "-in:file:silent protocols/jd3/pose_inputters/relaxed_decoys.out" );
		SilentFilePoseInputter inputter;
		PoseInputSources sources = inputter.pose_input_sources_from_command_line();
		TS_ASSERT_EQUALS( sources.size(), 4 );
		std::list< std::string > tags_list { "1amtFHA_0001", "1b0nFHB_0001", "1g1sFHC_0001", "1gjaFHA_0001" };
		utility::vector1< std::string > tags(4);
		std::copy( tags_list.begin(), tags_list.end(), tags.begin() );

		if ( sources.size() != 4 ) return;
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			TS_ASSERT_EQUALS( sources[ii]->origin(), SilentFilePoseInputter::keyname() );
			TS_ASSERT_EQUALS( sources[ii]->input_tag(), tags[ii] );
		}
	}

	void test_read_in_from_cmdline_without_pose_input_sources_call() {
		protocols_init_with_additional_options( "-in:file:silent protocols/jd3/pose_inputters/relaxed_decoys.out" );

		PoseInputSource input_source;
		input_source.origin( SilentFilePoseInputter::keyname() );
		input_source.input_tag( "1amtFHA_0001" );

		SilentFilePoseInputter inputter;

		core::pose::PoseOP pose = inputter.pose_from_input_source( input_source, basic::options::option, utility::tag::TagCOP() );

		TS_ASSERT_EQUALS( pose->annotated_sequence(), "P[PRO:NtermProteinFull]AAQVGLPVEQ[GLN:CtermProteinFull]" );
	}

};
