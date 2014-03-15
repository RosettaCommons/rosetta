// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/io/EmbedDefLoader.cxxtest.hh
///
/// @brief 		Test Suite for group of classes that load in an embedding definition resource
/// @details	Tests span file loader functionality using generic loader class
///
/// @author     Rebecca Alford

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/jd2/JD2ResourceManager.hh>

// Tested Classes
#include <core/membrane/io/EmbedDefLoader.hh>
#include <core/membrane/io/EmbedDefOptions.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

#include <protocols/loops/LoopsFileOptions.hh> //this is crazy - lol inspired by Milke :)

// Resource Manager Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/Job.hh>

// Basic
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

using namespace core::membrane::util;
using namespace core::membrane::io;
using namespace basic::resource_manager;
using namespace protocols::jd2;

/// @brief Test Class: Span File loader
class EmbedDefLoaderTests :  public CxxTest::TestSuite {
public: // public testing methods

	/// @brief SetUp - Runs before each test
	void setUp()
	{
        // Initialize
        protocols_init();
	}

	/// @brief tearDon - runs after each test
	void tearDown()
	{}

	/// Tests for IO Class ////

	/// @brief Test correct resource returned (spanning topology object)
	void test_returnType() {

        
        using namespace core::membrane::io;
        using namespace core::membrane::util;
        
        TS_TRACE("Testing create resource loader for embedding data");
        
        // Define opts and loader
        EmbedDefLoader loader;
        EmbedDefOptions opts;
        
        // Create an lstream
        std::string tag = "core/membrane/io/1afo_test.embed";
        std::istringstream lstream( tag );
        
        // Pull generic resource
        utility::pointer::ReferenceCountOP resource = loader.create_resource( opts, "core/membrane/io/1afo_test.embed", lstream );
        TS_ASSERT( resource );
        
        // Cast to an embed resource
        core::membrane::util::EmbedConfigInfoOP embed = dynamic_cast< EmbedConfigInfo * > ( resource () );
        TS_ASSERT( embed );

	}

	/// @brief Testing resource definition loading
	void test_resource_definition() {
        
        using namespace core::membrane::io;
        using namespace core::membrane::util;
        using namespace basic::resource_manager;
        using namespace protocols::jd2;
        
        // Create resource Manager Instances - Lazy Loading by Job for unit testing
        basic::resource_manager::ResourceManager * resource_manager( basic::resource_manager::ResourceManager::get_instance() );
		basic::resource_manager::LazyResourceManager * lazy_resource_manager( dynamic_cast< basic::resource_manager::LazyResourceManager * > ( resource_manager ));
        
        // Create a Job Stream from my membrane inputs
		std::istringstream a_jobstream(embed_def());
		JD2ResourceManagerJobInputter a_inputter;
		Jobs a_jobvector;
        
        // Try to fil the job
		try {
			a_inputter.fill_jobs_from_stream( a_jobstream, a_jobvector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "The following resource manager XML definitions file produced an exception:" << std::endl;
			std::cerr << "XML a:" << std::endl << embed_def() << std::endl;
			std::cerr << "Raised exception:" << std::endl << e.msg() << std::endl;
			lazy_resource_manager->show(std::cerr);
            
			TS_ASSERT( false );
		}
        
        // Load Mmebrane Embedding
        TS_TRACE("Loading per chain membrane embedding data");
        if (! lazy_resource_manager->has_resource_tag_by_job_tag("embedding", "membrane"))
        {
            throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
                                                    "-in:file:embedfile must be specified");
        }
        basic::resource_manager::ResourceOP embedding = lazy_resource_manager->get_resource_by_job_tag( "embedding", "membrane" );
        TS_ASSERT( embedding );

        // Checking that I returned the correct type
        core::membrane::util::EmbedConfigInfoOP embed = dynamic_cast< EmbedConfigInfo * > ( embedding () );
        TS_ASSERT( embed );
	}

private: // data
    
    /// @brief Job String
    std::string embed_def(){
		return
        "<JD2ResourceManagerJobInputter>\n"
        "  <ResourceLocators>\n"
        "    <FileSystemResourceLocator tag=1afo base_path=\"core/membrane/io/\" />\n"
        "  </ResourceLocators>\n"
        "  <Resources>\n"
        "    <PoseFromPDB tag=1afo_startstruct locator=1afo locatorID=\"1afo_test.pdb\" />\n"
        "    <EmbedDef tag=1afo_embed file=\"core/membrane/io/1afo_test.embed\" />\n"
        "  </Resources>\n"
        "  <Jobs>\n"
        "    <Job name=membrane>\n"
        "      <Data desc=\"startstruct\" resource_tag=1afo_startstruct />\n"
        "      <Data desc=\"embedding\" resource_tag=1afo_embed />\n"
        "    </Job>\n"
        "  </Jobs>\n"
        "</JD2ResourceManagerJobInputter>\n";
	}
};
