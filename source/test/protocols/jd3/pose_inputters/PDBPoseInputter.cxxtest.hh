// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/JD2ResourceManager.cxxtest.hh
/// @brief  test suite for protocols::jd2::JD2ResourceManager and protocols::resource_manager::LazyResourceManager
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <protocols/jd3/pose_inputters/PDBPoseInputter.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.hh>

#include <basic/options/util.hh>

#include <utility/options/keys/OptionKeyList.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>

template <typename T>
typename T::mapped_type get(T const& map, typename T::key_type const& key)
{
	typename T::const_iterator iter(map.find(key));
	return iter != map.end() ? iter->second : typename T::mapped_type();
}

using protocols::jd3::pose_inputters::PoseInputSources;
using protocols::jd3::pose_inputters::PoseInputSourceOP;
using protocols::jd3::pose_inputters::PDBPoseInputter;

class JD3PDBPoseInputterTests : public CxxTest::TestSuite {
public:

	void setUp() {
	}

	utility::vector1< std::string >
	four_pdbs_list_contents() {
		std::list< std::string > pdbs = { "1abc.pdb", "2def.pdb", "3ghi.pdb", "4jkl.pdb" };
		utility::vector1< std::string > pdbs_vect( 4 );
		std::copy( pdbs.begin(), pdbs.end(), pdbs_vect.begin() );
		return pdbs_vect;
	}


	void test_read_s_flag_one_pdb() {
		protocols_init_with_additional_options( "-s /home/andrew/1ubq.pdb" );
		PDBPoseInputter inputter;
		PoseInputSources sources = inputter.pose_input_sources_from_command_line();
		TS_ASSERT_EQUALS( sources.size(), 1 );
		if ( sources.size() != 1 ) return;
		TS_ASSERT_EQUALS( sources[1]->input_tag(), std::string("1ubq") );
		TS_ASSERT_EQUALS( sources[1]->origin(), "PDB" );
		TS_ASSERT( sources[1]->string_string_map().find( std::string("filename") ) != sources[1]->string_string_map().end() );
		TS_ASSERT_EQUALS( get( sources[1]->string_string_map(), std::string("filename")), "/home/andrew/1ubq.pdb" );
	}

	void test_read_s_flag_three_pdbs() {
		protocols_init_with_additional_options( "-s /home/andrew/1ubq.pdb /net/scr/1ubi.pdb /www/whatever/1l2y.pdb" );
		PDBPoseInputter inputter;
		PoseInputSources sources = inputter.pose_input_sources_from_command_line();
		TS_ASSERT_EQUALS( sources.size(), 3 );
		if ( sources.size() != 3 ) return;

		TS_ASSERT_EQUALS( sources[1]->input_tag(), std::string("1ubq") );
		TS_ASSERT_EQUALS( sources[1]->origin(), "PDB" );
		TS_ASSERT( sources[1]->string_string_map().find( std::string("filename") ) != sources[1]->string_string_map().end() );
		TS_ASSERT_EQUALS( get( sources[1]->string_string_map(), std::string("filename")), "/home/andrew/1ubq.pdb" );

		TS_ASSERT_EQUALS( sources[2]->input_tag(), std::string("1ubi") );
		TS_ASSERT_EQUALS( sources[2]->origin(), "PDB" );
		TS_ASSERT( sources[2]->string_string_map().find( std::string("filename") ) != sources[2]->string_string_map().end() );
		TS_ASSERT_EQUALS( get( sources[2]->string_string_map(), std::string("filename")), "/net/scr/1ubi.pdb" );

		TS_ASSERT_EQUALS( sources[3]->input_tag(), std::string("1l2y") );
		TS_ASSERT_EQUALS( sources[3]->origin(), "PDB" );
		TS_ASSERT( sources[3]->string_string_map().find( std::string("filename") ) != sources[3]->string_string_map().end() );
		TS_ASSERT_EQUALS( get( sources[3]->string_string_map(), std::string("filename")), "/www/whatever/1l2y.pdb" );
	}

	void test_read_l_flag_four_pdbs() {
		protocols_init_with_additional_options( "-l protocols/jd3/pose_inputters/four_pdbs.list" );
		PDBPoseInputter inputter;
		PoseInputSources sources = inputter.pose_input_sources_from_command_line();
		TS_ASSERT_EQUALS( sources.size(), 4 );
		if ( sources.size() != 4 ) return;

		utility::vector1< std::string > pdbs = four_pdbs_list_contents();
		for ( core::Size ii = 1; ii <= sources.size(); ++ii ) {
			TS_ASSERT_EQUALS( pdbs[ ii ].substr( 0, 4 ), sources[ ii ]->input_tag() );
			TS_ASSERT_EQUALS( get( sources[ii]->string_string_map(), std::string("filename")), pdbs[ii] );
		}
	}

	void test_read_input_sources_from_listfile_attribute_in_tag() {
		using namespace utility::tag;

		protocols_init_with_additional_options( "" );
		utility::options::OptionKeyList opts_list;
		PDBPoseInputter::list_options_read( opts_list );
		utility::options::OptionCollectionOP opts = basic::options::read_subset_of_global_option_collection( opts_list );

		PDBPoseInputter inputter;
		TagCOP input_tag = Tag::create( "<PDB listfile=\"protocols/jd3/pose_inputters/four_pdbs.list\"/>" );

		PoseInputSources sources = inputter.pose_input_sources_from_tag( *opts, input_tag );

		TS_ASSERT_EQUALS( sources.size(), 4 );
		if ( sources.size() != 4 ) return;

		utility::vector1< std::string > pdbs = four_pdbs_list_contents();
		for ( core::Size ii = 1; ii <= sources.size(); ++ii ) {
			TS_ASSERT_EQUALS( pdbs[ ii ].substr( 0, 4 ), sources[ ii ]->input_tag() );
			TS_ASSERT_EQUALS( get( sources[ii]->string_string_map(), std::string("filename")), pdbs[ii] );
		}
	}

	void test_read_input_sources_from_listfile_attribute_in_tag_w_path() {
		using namespace utility::tag;

		protocols_init_with_additional_options( "" );
		utility::options::OptionKeyList opts_list;
		PDBPoseInputter::list_options_read( opts_list );
		utility::options::OptionCollectionOP opts = basic::options::read_subset_of_global_option_collection( opts_list );

		PDBPoseInputter inputter;
		TagCOP input_tag = Tag::create( "<PDB listfile=\"protocols/jd3/pose_inputters/four_pdbs.list\" path=\"protocols/jd3/pose_inputters\"/>" );

		PoseInputSources sources = inputter.pose_input_sources_from_tag( *opts, input_tag );

		TS_ASSERT_EQUALS( sources.size(), 4 );
		if ( sources.size() != 4 ) return;

		utility::vector1< std::string > pdbs = four_pdbs_list_contents();
		for ( core::Size ii = 1; ii <= sources.size(); ++ii ) {
			TS_ASSERT_EQUALS( pdbs[ ii ].substr( 0, 4 ), sources[ ii ]->input_tag() );
			TS_ASSERT_EQUALS( get( sources[ii]->string_string_map(), std::string("filename")), "protocols/jd3/pose_inputters/" + pdbs[ii] );
		}
	}

};
