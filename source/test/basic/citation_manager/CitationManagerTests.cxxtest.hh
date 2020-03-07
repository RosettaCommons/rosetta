// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  basic/citation_manager/CitationManagerTests.cxxtest.hh
/// @brief  Unit tests for the citation manager, a class for tracking what modules have been
/// used during a Rosetta run and what papers or authors should be cited.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/Citation.hh>

static basic::Tracer TR("CitationManagerTests");


class CitationManagerTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();
	}

	void tearDown(){ }

	/// @brief Test that we've loaded citations from the database correctly.
	void test_database_load() {
		using namespace basic::citation_manager;
		CitationCOP pnaspaper( CitationManager::get_instance()->get_citation_by_doi( "10.1073/pnas.1710695114" ) );
		TS_ASSERT_EQUALS( pnaspaper->article_title(), "De novo design of covalently constrained mesosize protein scaffolds with unique tertiary structures" );
		TS_ASSERT_EQUALS( pnaspaper->journal_title(), "Proc Natl Acad Sci USA" );
		TS_ASSERT_EQUALS( pnaspaper->year(), 2017 );
		TS_ASSERT_EQUALS( pnaspaper->volume_issue_pages(), "114(41):10852–10857" );
		TS_ASSERT_EQUALS( pnaspaper->senior_authors().size(), 1 );
		TS_ASSERT_EQUALS( pnaspaper->senior_authors()[1].given_names(), "William" );
		TS_ASSERT_EQUALS( pnaspaper->senior_authors()[1].surname(), "DeGrado" );
		TS_ASSERT_EQUALS( pnaspaper->senior_authors()[1].initials(), "WF");
		TS_ASSERT_EQUALS( pnaspaper->co_authors().size(), 6 );
		TS_ASSERT_EQUALS( pnaspaper->co_authors()[5].given_names(), "Daniel-Adriano" );
		TS_ASSERT_EQUALS( pnaspaper->co_authors()[5].surname(), "Silva" );
		TS_ASSERT_EQUALS( pnaspaper->co_authors()[5].initials(), "D-A" );
		TS_ASSERT_EQUALS( pnaspaper->primary_authors().size(), 3 );
		TS_ASSERT_EQUALS( pnaspaper->primary_authors()[2].surname(), "Wu" );
		TS_ASSERT_EQUALS( pnaspaper->primary_authors()[1].given_names(), "Bobo" );
		TS_ASSERT_EQUALS( pnaspaper->primary_authors()[3].initials(), "VK" );

		std::stringstream outstream;
		TS_ASSERT_THROWS_NOTHING( pnaspaper->get_formatted_citation( outstream ) );
		TR << outstream.str() << std::endl;
	}


};
