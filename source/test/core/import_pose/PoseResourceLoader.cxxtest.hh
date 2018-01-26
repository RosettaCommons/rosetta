// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/import_pose/PoseResourceLoader.cxxtest.hh
/// @brief test suite for core::import::PoseResourceLoader
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/import_pose/PoseResourceLoader.hh>
#include <core/import_pose/PoseResource.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Package headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <sstream>

using namespace basic::resource_manager;

class PoseFromPDBLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}

	// @brief test default options and default locator
	void test_PoseFromPDBLoader() {

		using namespace core::import_pose;

		std::string test_pdb_filename = "core/io/test_in_idealized.pdb";
		std::string contents = utility::file_contents( test_pdb_filename );
		std::istringstream iss1( contents );
		std::istringstream iss2( contents );
		utility::tag::TagCOP tag( utility::tag::Tag::create( "<Pose><PDB/></Pose>" ));

		ResourceManager rm;
		PoseResourceLoader loader;

		ResourceCOP my_resource = loader.create_resource( rm,  tag, "unit test", iss1 );
		PoseResourceCOP pose_resource = utility::pointer::dynamic_pointer_cast< PoseResource const > ( my_resource );
		core::pose::PoseCOP idealized_input_pose = pose_resource->pose_deep_copy();

		TS_ASSERT( idealized_input_pose ); // make sure we got back the right resource type

		// Crearte a pose in the old fashioned way and make sure the contents of this pose and
		// the pose from the Loader are identical.
		core::pose::Pose reference_pose;
		ImportPoseOptions opts;
		core::import_pose::pose_from_pdb_stream( reference_pose, iss2, "unit test", opts );

		TS_ASSERT(
			core::pose::compare_binary_protein_silent_struct(
			*idealized_input_pose, reference_pose));

	}

};
