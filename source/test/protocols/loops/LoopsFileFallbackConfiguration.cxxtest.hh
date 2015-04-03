// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsFileOptions.cxxtest.hh
/// @brief test suite for protocols/loops/LoopsFileOptions.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/loops/LoopsFileFallbackConfiguration.hh>
#include <basic/resource_manager/FallbackConfigurationFactory.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/types.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Numeric headers
//#include <numeric/kdtree/WrappedPrimitive.hh>

// C++ headers
#include <string>


class LoopsFileFallbackConfigurationTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	void test_use_FallbackConfigFactory_to_get_loops_file() {
		using namespace protocols::loops;
		using namespace basic::resource_manager;

		bool better_be_true = FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( "LoopsFile" );
		TS_ASSERT( better_be_true );

		bool better_be_false = FallbackConfigurationFactory::get_instance()->has_fallback_for_resource( "nonsense!" );
		TS_ASSERT( ! better_be_false );

		std::string loop_filename = "this_could_be_anything";

		try {
			FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( "LoopsFile" );
			std::string loopfilename = fallback->get_locator_id("LoopsFile" );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "The fallback LoopsFile resource option has no loops files associated with it! Was the option omitted from the command line?";
			TS_ASSERT( expected_error_message == e.msg() );
		}

		basic::options::option[ basic::options::OptionKeys::loops::loop_file ].def( loop_filename );
		FallbackConfigurationOP fallback = FallbackConfigurationFactory::get_instance()->create_fallback_configuration( "LoopsFile" );

		LoopsFileFallbackConfigurationOP loops_fallback = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileFallbackConfiguration > ( fallback );
		TS_ASSERT( loops_fallback ); // make sure the cast worked
		TS_ASSERT( loop_filename == loops_fallback->get_locator_id( "LoopsFile" ) );
		TS_ASSERT( "LoopsFile" == loops_fallback->get_resource_loader( "LoopsFile" ) );
		TS_ASSERT( 0 == loops_fallback->get_resource_options( "LoopsFile" ) );

	}

};
