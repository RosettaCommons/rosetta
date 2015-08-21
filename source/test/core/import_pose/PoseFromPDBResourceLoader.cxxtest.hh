// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/import_pose/PoseFromPDBResourceLoader.cxxtest.hh
/// @brief test suite for core::import::PoseFromPDBLoader
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Package headers
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/ResourceManagerFactory.hh>
#include <basic/resource_manager/types.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>

// Numberic headers

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

		std::string test_pdb_filename = "core/io/test_in_idealized.pdb";

		ResourceTag rTag = "an_ideal_pose_for_testing";
		ResourceDescription rDesc = "idealized_input_pose";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = test_pdb_filename;
		my_config.loader_type = "PoseFromPDB";
		my_config.resource_options_tag = "";

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type
		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
		core::pose::PoseCOP idealized_input_pose = utility::pointer::dynamic_pointer_cast< core::pose::Pose const > ( my_resource );

		TS_ASSERT( idealized_input_pose ); // make sure we got back the right resource type

		// Crearte a pose in the old fashioned way and make sure the contents of this pose and
		// the pose from the Loader are identical.
		core::pose::Pose reference_pose;
		core::import_pose::pose_from_pdb( reference_pose, test_pdb_filename );

		TS_ASSERT(
			core::pose::compare_binary_protein_silent_struct(
			*idealized_input_pose, reference_pose));

	}

};
