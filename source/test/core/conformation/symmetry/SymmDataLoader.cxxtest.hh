// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/conformation/symmetry/SymmDataLoader.cxxtest.hh
/// @brief test suite for core::conformation::symmetry::SymmDataLoader
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmDataLoader.hh>

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

class SymmDataLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
	}

	void test_SymmDataLoader(){
	  std::string symm_def1_filename("core/conformation/symmetry/symm_def1.dat");
		do_test_SymmDataLoader(symm_def1_filename);
	}

	// @brief test default options and default locator
	void do_test_SymmDataLoader(
		std::string const & symm_data_filename
) {
		// Prepare the resource manager so we can ask it for the symmetry
		// definition
		ResourceTag rTag = "symm_def1";
		ResourceDescription rDesc = "SymmData";
		JobTag jTag = "1";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = symm_data_filename;
		my_config.loader_type = "SymmData";
		my_config.resource_options_tag = "";

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager(
			dynamic_cast< LazyResourceManager * > ( resource_manager ));

		// make sure we got back the right resource_manager type
		TS_ASSERT( lazy_resource_manager );
		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		ResourceOP my_resource(
			lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag ));
		core::conformation::symmetry::SymmDataCOP symm_data(
			utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmData const > ( my_resource ));

		// make sure we got back the right resource type
		TS_ASSERT( symm_data );

		// Setup the alternate symmetry data to compare against
		core::conformation::symmetry::SymmDataOP symm_data_alt( new core::conformation::symmetry::SymmData() );
		symm_data_alt->read_symmetry_data_from_file(symm_data_filename);

		TS_ASSERT_EQUALS(*symm_data, *symm_data_alt);
	}

};
