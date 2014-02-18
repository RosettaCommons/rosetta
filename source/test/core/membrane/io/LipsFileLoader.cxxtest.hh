// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/SpanFileIO.cxxtest.hh
///
/// @brief 		Test Suite for izstream reader class that reads in OCTOPUS file data
/// @details
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Tested Classes
#include <core/membrane/io/LipoFileLoader.hh>

#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/membrane/io/LipoFileOptions.hh>
#include <protocols/loops/LoopsFileOptions.hh>

// Resource Manager Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// Package Headers
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class: Lipo File loader
class LipoFileLoaderTest : public CxxTest::TestSuite {
    
public: // test methods
    
	/// @brief SetUp - Runs before each test
	void setUp()
	{
        
        using namespace core::membrane::io;
        using namespace protocols::loops;
        
        // Core init
        core_init();
	}
    
	/// @brief tearDon - runs after each test
	void tearDown()
	{}
    
	/// Tests for IO Class ////
    /// @brief Test correct response to an empty string locator
    void test_emptyLocator()
    {
        // Empty Lstream
        std::istringstream lstream("unit_test");
        
        TS_TRACE("Testing empty string locator - invalid case...");
        TS_ASSERT_THROWS_ANYTHING( loader_.create_resource( opts1_, "", lstream ) );
        
    }
    
    /// @brief Test correct resource returned (spanning topology object)
    void test_returnType()
    {
        
        using namespace core::membrane::util;
        using namespace core::membrane::properties;
        
        // Empty Lstream
        std::istringstream lstream("unit_test");
        
        // Load the data
        utility::pointer::ReferenceCountOP lips_exp = loader_.create_resource( opts1_, "core/membrane/io/1afo_test.lips", lstream );

        // Ensure resources were returned
        TS_ASSERT( lips_exp );
        
        // Cast to type
        core::membrane::properties::LipidAccInfoOP lpptr1 = dynamic_cast< core::membrane::properties::LipidAccInfo * > ( lips_exp () );
        
        // Check non null (typechecking)
        TS_ASSERT( lpptr1 );
    }
    
    /// @brief Testing resource definition
    void test_resource_definition()
    {
        
        using namespace core::membrane::io;
        using namespace core::membrane::util;
        using namespace core::membrane::properties;
        using namespace basic::resource_manager;
        using namespace protocols::jd2;
        
        // Create resource Manager Instances - Lazy Loading by Job for unit testing
        basic::resource_manager::ResourceManager * resource_manager( basic::resource_manager::ResourceManager::get_instance() );
		basic::resource_manager::LazyResourceManager * lazy_resource_manager( dynamic_cast< basic::resource_manager::LazyResourceManager * > ( resource_manager ));
        
        // Create a Job Stream from my membrane inputs
		std::istringstream a_jobstream(lips_def());
		JD2ResourceManagerJobInputter a_inputter;
		Jobs a_jobvector;
        
        // Try to fil the job
		try {
			a_inputter.fill_jobs_from_stream( a_jobstream, a_jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML a:" << std::endl << lips_def() << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);
            
			TS_ASSERT( false );
		}
        
        // Load Membrane Lips File
        TS_TRACE("Loading a membrane lipids obj from the resource manager");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("lipids", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:lipofile must be specified");
        }
        basic::resource_manager::ResourceOP lipids = lazy_resource_manager->get_resource_by_job_tag( "lipids", "membrane" );
        TS_ASSERT( lipids );
        
        // Cast to type
        core::membrane::properties::LipidAccInfoOP lpptr1 = dynamic_cast< LipidAccInfo * > ( lipids () );
        
        // Check non null (typechecking)
        TS_ASSERT( lpptr1 );
    }

private: // resource definition
    
    /// @brief Resource definition for lipid defs
    std::string lips_def() {
        return
        "<JD2ResourceManagerJobInputter>\n"
        "  <ResourceLocators>\n"
        "    <FileSystemResourceLocator tag=1afo base_path=\"core/membrane/io/\" />\n"
        "  </ResourceLocators>\n"
        "  <Resources>\n"
        "    <PoseFromPDB tag=1afo_startstruct locator=1afo locatorID=\"1afo_test.pdb\" />\n"
        "    <LipoFile tag=1afo_lips file=\"core/membrane/io/1afo_test.lips\" />\n"
        "  </Resources>\n"
        "  <Jobs>\n"
        "    <Job name=membrane>\n"
        "      <Data desc=\"startstruct\" resource_tag=1afo_startstruct />\n"
        "      <Data desc=\"lipids\" resource_tag=1afo_lips />\n"
        "    </Job>\n"
        "  </Jobs>\n"
        "</JD2ResourceManagerJobInputter>\n";
    }
    
private: // data
    
    // Set up a permanent loader class (gets reinitialized per test)
    core::membrane::io::LipoFileLoader loader_;
    
    // Storing opts pointers
    core::membrane::io::LipoFileOptions opts1_;
    protocols::loops::LoopsFileOptions opts2_; // this is hysterical
    
};
