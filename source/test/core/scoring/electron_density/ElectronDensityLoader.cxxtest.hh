// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/import_pose/ElectronDensityResourceLoader.cxxtest.hh
/// @brief test suite for core::import::ElectronDensityLoader
/// @author Matthew O'Meara mattjomeara@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Unit headers
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/ElectronDensityLoader.hh>

// Project headers
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/ResourceManagerFactory.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>
#include <basic/resource_manager/types.hh>

// Package headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/EXCN_Base.hh>

// C++ headers
#include <sstream>

static basic::Tracer TR("test.scoring.electron_density.ElectronDensityLoader");

using namespace basic::resource_manager;

class ElectronDensityLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}


	// @brief test default options and default locator
	void test_ElectronDensityLoader() {

		try{
			do_test_ElectronDensityLoader();

		} catch( utility::excn::EXCN_Base& excn ) {
			TR
				<< "ERROR: Exception caught by ElectronDensityLoaderTests:"
				<< excn << std::endl;
		}
	}
	void do_test_ElectronDensityLoader() {
		
		std::string test_electron_density_filename =
			"core/scoring/electron_density/1onc_8A.mrc";

		ResourceTag rTag = "electron_density_tag";
		ResourceDescription rDesc = "electron_density_description";
		JobTag jTag = "job_tag";

		ResourceLocatorOP binary_file_locator( new locator::FileSystemResourceLocator(
				std::ios_base::binary | std::ios_base::in ) );

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "binary_file";
		my_config.locator_id = test_electron_density_filename;
		my_config.loader_type = "ElectronDensityLoader";
		my_config.resource_options_tag = "";

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager =
			dynamic_cast< LazyResourceManager * > ( resource_manager );

		// make sure we got back the right resource_manager type
		TS_ASSERT( lazy_resource_manager );

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		lazy_resource_manager->add_resource_locator(
			my_config.locator_tag, binary_file_locator);

		ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag(
			rDesc, jTag );
		core::scoring::electron_density::ElectronDensityOP electron_density(
			utility::pointer::dynamic_pointer_cast< core::scoring::electron_density::ElectronDensity >(my_resource));

		// make sure we got back the right resource type
		TS_ASSERT( electron_density );

		electron_density->writeMRC("1onc_8A_out.mrc");

		//TODO figure out how to test if the electron density file has been loaded correctly.
	}

};
