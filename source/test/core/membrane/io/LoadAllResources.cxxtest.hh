// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/LoadAllResourcesTest.cxxtest.hh
///
/// @brief 		Test Suite for loading bottom level memrbane resources
/// @details    This test suite should act as a top level test that membrane resources can eb grabbed form
///             an instance of the resource manager when associated with a given job.
///
/// @note       Last Modified: 12/20/13
/// @author 	Rebecca Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>
#include <core/pose/Pose.hh>

// Package Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

/// @brief Test Class Load All Membrane Resources
/// @details Load all membrane resources - example code for unit tests
class LoadAllResourcesTest :  public CxxTest::TestSuite {
public:
    
    /// @brief Set Up
    void setUp() {
        protocols_init();
    }
    
    /// @brief Tear down
    void tearDown()
    {}

    /// @brief Test Loading Membrane Resources
    void test_load_membrane_resources() {
        
        TS_TRACE("Testing loading membrane resources");
        
        using namespace protocols::jd2;
        using namespace core::membrane::util;

        // Create resource Manager Instances - Lazy Loading by Job for unit testing
        basic::resource_manager::ResourceManager * resource_manager( basic::resource_manager::ResourceManager::get_instance() );
		basic::resource_manager::LazyResourceManager * lazy_resource_manager( dynamic_cast< basic::resource_manager::LazyResourceManager * > ( resource_manager ));
        
        // Create a Job Stream from my membrane inputs
		std::istringstream a_jobstream(membrane_input());
		JD2ResourceManagerJobInputter a_inputter;
		Jobs a_jobvector;
        
        // Try to fil the job
		try {
			a_inputter.fill_jobs_from_stream( a_jobstream, a_jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML a:" << std::endl << membrane_input() << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);
            
			TS_ASSERT( false );
		}
        
        TS_TRACE("Showing the current resources in my JD2 resource manager: ");
        lazy_resource_manager->show(std::cout);
        
        
        // Load Pose
        TS_TRACE("Loading a pose from the resource manager");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("startstruct", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:s must be specified");
        }
        basic::resource_manager::ResourceOP pose = lazy_resource_manager->get_resource_by_job_tag( "startstruct", "membrane" );
        TS_ASSERT( pose );
        
        // Load Membrane topology
        TS_TRACE("Loading a membrane topology obj from the resource manager");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("topology", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:spanfile must be specified");
        }
        basic::resource_manager::ResourceOP topology = lazy_resource_manager->get_resource_by_job_tag( "topology", "membrane" );
        TS_ASSERT( topology );
        
        // Load Mmebrane Embedding
        TS_TRACE("Loading per chain membrane embedding data");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("embedding", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:embedfile must be specified");
        }
        basic::resource_manager::ResourceOP embedding = lazy_resource_manager->get_resource_by_job_tag( "embedding", "membrane" );
        TS_ASSERT( embedding );
        
        // Loading Membrane Embedding Options
        TS_TRACE("loading membrane embedding parameters");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("params", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception("Either a resource definition file or the command line option -in:membrane:params "\
                                                    "must be specified");
        }
        basic::resource_manager::ResourceOP params = lazy_resource_manager->get_resource_by_job_tag( "params", "membrane" );
        TS_ASSERT( params );
        
        // Load Membrane Lips File
        TS_TRACE("Loading a membrane lipids obj from the resource manager");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("lipids", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:lipofile must be specified");
        }
        basic::resource_manager::ResourceOP lipids = lazy_resource_manager->get_resource_by_job_tag( "lipids", "membrane" );
        TS_ASSERT( lipids );
        
        
        TS_TRACE("Test successful!");
    }
    
    /// @brief Job String
    std::string membrane_input(){
		return
        "<JD2ResourceManagerJobInputter>\n"
        "  <ResourceLocators>\n"
        "    <FileSystemResourceLocator tag=1afo base_path=\"core/membrane/io/\" />\n"
        "  </ResourceLocators>\n"
        "  <ResourceOptions>\n"
        "    <EmbedSearchParamsOptions tag=embed_options normal_search=1 />\n"
        "  </ResourceOptions>\n"
        "  <Resources>\n"
        "    <PoseFromPDB tag=1afo_startstruct locator=1afo locatorID=\"1afo_test.pdb\" />\n"
        "    <SpanFile tag=1afo_topology file=\"core/membrane/io/1afo_test.span\" />\n"
        "    <EmbedDef tag=1afo_embed file=\"core/membrane/io/1afo_test.embed\" />\n"
        "    <EmbedSearchParams tag=1afo_params file=\"core/membrane/io/1afo_test.embed\" />\n"
        "    <LipoFile tag=1afo_lips file=\"core/membrane/io/1afo_test.lips\" />\n"
        "  </Resources>\n"
        "  <Jobs>\n"
        "    <Job name=membrane>\n"
        "      <Data desc=\"startstruct\" resource_tag=1afo_startstruct />\n"
        "      <Data desc=\"topology\" resource_tag=1afo_topology />\n"
        "      <Data desc=\"embedding\" resource_tag=1afo_embed />\n"
        "      <Data desc=\"params\" resource_tag=1afo_params />\n"
        "      <Data desc=\"lipids\" resource_tag=1afo_lips />\n"
        "    </Job>\n"
        "  </Jobs>\n"
        "</JD2ResourceManagerJobInputter>\n";
	}
};
