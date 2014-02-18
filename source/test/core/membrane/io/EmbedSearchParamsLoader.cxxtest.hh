// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/EmbedSearchParamsLoader.cxxtest.hh
///
/// @brief 		Test Suite for loading an embedding search params loader
/// @details
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Tested Classes
#include <core/membrane/io/EmbedSearchParamsLoader.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <protocols/loops/LoopsFileOptions.hh>
#include <core/membrane/io/EmbedSearchParamsOptions.hh>

// Resource Manager Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace core::membrane::util;
using namespace core::membrane::io;

/// @brief Test Class: Span File loader
class EmbedSearchParamsLoaderTests :  public CxxTest::TestSuite {

public: // public testing methods

	/// @brief SetUp - Runs before each test
	void setUp()
	{
        // setup
        core_init();
	}

	/// @brief tearDon - runs after each test
	void tearDown()
	{}

	/// Tests for IO Class ////

	/// @brief Test correct resource returned (spanning topology object)
	void test_returnType()
	{
		TS_TRACE("Testing create_resource method");
        
		// Declare a new loader and otps
		core::membrane::io::EmbedSearchParamsLoader loader;
		core::membrane::io::EmbedSearchParamsOptionsOP esp_opts_ = new core::membrane::io::EmbedSearchParamsOptions();
        
		// Create a ut lstream
		std::istringstream lstream;
        
		// Initializing esp opts
		esp_opts_->set_normal_search(true);
		esp_opts_->set_normal_start_angle(10);
		esp_opts_->set_normal_delta_angle(10);
		esp_opts_->set_normal_max_angle(10);
		esp_opts_->set_center_search(true);
		esp_opts_->set_center_max_delta(10);
		esp_opts_->set_normal_mag(20);
		esp_opts_->set_center_mag(30);
		esp_opts_->set_normal_cycles(10);
		esp_opts_->set_penalties(true);
		esp_opts_->set_no_interpolate_mpair(true);
        
		// Create opts ref obj
		core::membrane::io::EmbedSearchParamsOptions const & opts ( *esp_opts_ );
        
		utility::pointer::ReferenceCountOP resource = loader.create_resource( opts, "", lstream );
	}

	/// @brief Testing resource definition from rm
	void test_resource_definition()
	{

        using namespace core::membrane::io;
        using namespace core::membrane::util;
        using namespace basic::resource_manager;
        using namespace protocols::jd2;
        
        // Create resource Manager Instances - Lazy Loading by Job for unit testing
        basic::resource_manager::ResourceManager * resource_manager( basic::resource_manager::ResourceManager::get_instance() );
		basic::resource_manager::LazyResourceManager * lazy_resource_manager( dynamic_cast< basic::resource_manager::LazyResourceManager * > ( resource_manager ));
        
        // Create a Job Stream from my membrane inputs
		std::istringstream a_jobstream(params_def());
		JD2ResourceManagerJobInputter a_inputter;
		Jobs a_jobvector;
        
        // Try to fil the job
		try {
			a_inputter.fill_jobs_from_stream( a_jobstream, a_jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML a:" << std::endl << params_def() << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);
            
			TS_ASSERT( false );
		}
        
        // Loading Membrane Embedding Options
        TS_TRACE("loading membrane embedding parameters");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("params", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception("Either a resource definition file or the command line option -in:membrane:params "\
                                                    "must be specified");
        }
        basic::resource_manager::ResourceOP params = lazy_resource_manager->get_resource_by_job_tag( "params", "membrane" );
        TS_ASSERT( params );
        
	}

private: // resource definition
    
    /// @brief Job Definition input (resource definition)
    std::string params_def(){
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
        "    <EmbedSearchParams tag=1afo_params locator=1afo locatorID=\"1afo_test.embed\" />\n"
        "  </Resources>\n"
        "  <Jobs>\n"
        "    <Job name=membrane>\n"
        "      <Data desc=\"startstruct\" resource_tag=1afo_startstruct />\n"
        "      <Data desc=\"params\" resource_tag=1afo_params />\n"
        "    </Job>\n"
        "  </Jobs>\n"
        "</JD2ResourceManagerJobInputter>\n";
	}

};
