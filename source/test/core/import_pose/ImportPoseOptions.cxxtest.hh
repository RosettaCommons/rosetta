// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/import_pose/ImportPoseOptions.cxxtest.hh
/// @brief  Test suite for ImportPoseOptions's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/import_pose/import_pose_options.hh>

// Project headers
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/backtrace.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>

// basic headers
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh> // for a handful of options that should be read
#include <basic/options/keys/packing.OptionKeys.gen.hh> // for a handful of options that should not be read

// C++ headers
#include <algorithm>

class ImportPoseOptionsTests : public CxxTest::TestSuite {
public:
	typedef utility::keys::VariantKey< utility::options::OptionKey > VariantOptionKey;

	/// @brief Make sure that ImportPoseOptions::list_options_read lists every one of
	/// the options that the ImportPoseOptions object reads in its init_from_options
	/// method.  Bad things happen if these two functions get out of sync.
	void test_options_and_options_list_in_sync()
	{

		core_init();
		utility::options::OptionKeyList import_pose_options;
		core::import_pose::ImportPoseOptions::list_options_read( import_pose_options );
		// cursory check of some of the options known to be read by the ImportPoseOptions
		TS_ASSERT_EQUALS( std::count( import_pose_options.begin(), import_pose_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::centroid )), 1 );
		TS_ASSERT_EQUALS( std::count( import_pose_options.begin(), import_pose_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::fullatom )), 1 );
		TS_ASSERT_EQUALS( std::count( import_pose_options.begin(), import_pose_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::residue_type_set )), 1 );

		// cursory check that not all options are in here, because that would be weird
		TS_ASSERT_EQUALS( std::count( import_pose_options.begin(), import_pose_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::ex1::ex1 )), 0 );

		utility::options::OptionCollectionCOP import_pose_option_collection =
			basic::options::read_subset_of_global_option_collection( import_pose_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure(); // just in case
			core::import_pose::ImportPoseOptions import_pose( *import_pose_option_collection );
		} catch ( ... ) {
			TS_ASSERT( false ); // we screwed the pooch
		}
	}

	void test_options_actually_reads_option_collection()
	{
		core_init();
		utility::options::OptionKeyList import_pose_options;
		core::import_pose::ImportPoseOptions::list_options_read( import_pose_options );

		// Now drop one of the options from the list, and make sure that when we construct an
		// ImportPoseOptions from the option collection that an exception gets thrown.
		import_pose_options.pop_front();

		utility::options::OptionCollectionCOP import_pose_option_collection =
			basic::options::read_subset_of_global_option_collection( import_pose_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure();
			core::import_pose::ImportPoseOptions import_pose( *import_pose_option_collection );
			TS_ASSERT( false ); // we screwed the pooch
		} catch ( ... ) {
			// good -- if we don't list an option that we're going to read, then
			// an exception will be thrown / an assertion failure will get triggered
			TS_ASSERT( true );
		}

	}

};
