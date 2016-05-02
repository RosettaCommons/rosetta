// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief  tests for the Antibody Design options classes and parsers
/// @author Jared Adolf-Bryfogle


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/database/CDRSetOptionsParser.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

// Utility headers
#include <utility/backtrace.hh>
#include <utility/vector1.hh>

#define BFE BOOST_FOREACH
using namespace protocols::antibody;
using namespace protocols::antibody::design;
using namespace protocols::antibody::clusters;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.AntibodyDatabase");

class AntibodyDatabase: public CxxTest::TestSuite {
	core::pose::Pose pose;
	AntibodyInfoOP ab_info;


public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_file(pose, "protocols/antibody/aho_with_antigen.pdb", core::import_pose::PDB_file); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );

	}

	void test_options_cdr_set(){


		CDRSetOptionsParserOP parser( new CDRSetOptionsParser() );
		utility::vector1< CDRSetOptionsOP > options = parser->parse_options("protocols/antibody/database/cdr_set_instructions.txt");

		utility::vector1< std::string > include_species;
		include_species.push_back( "Hu" );

		utility::vector1< std::string > exclude_species;
		exclude_species.push_back( "Mo" );

		utility::vector1< std::string > include_germlines;
		include_germlines.push_back( "IGKV1" );

		utility::vector1< std::string > exclude_germlines;
		exclude_germlines.push_back( "IGKV3" );

		utility::vector1< CDRClusterEnum > include_l1_clusters;
		include_l1_clusters.push_back( L1_11_1 );

		utility::vector1< CDRClusterEnum > exclude_l1_clusters;
		exclude_l1_clusters.push_back( L1_11_1 );

		utility::vector1< std::string > include_pdbs;
		include_pdbs.push_back( "1HIL" );
		include_pdbs.push_back( "2J88" );

		utility::vector1< std::string > exclude_pdbs;
		exclude_pdbs.push_back( "1HIL" );

		utility::vector1< bool > length_types( 3, false );
		length_types[ 1 ] = true;

		for ( core::Size i = 1; i <= 6; ++i ) {
			CDRSetOptionsOP cdr_options = options[ i ];
			TS_ASSERT_EQUALS( cdr_options->length_type(), length_types );
			TS_ASSERT_EQUALS( cdr_options->max_length(), 30 );
			TS_ASSERT_EQUALS( cdr_options->min_length(), 1 );

			TS_ASSERT_EQUALS( cdr_options->include_only_center_clusters(), true);
			TS_ASSERT_EQUALS( cdr_options->include_only_current_cluster(), true);

			TS_ASSERT_EQUALS( cdr_options->include_only_species(), include_species );
			TS_ASSERT_EQUALS( cdr_options->exclude_species(), exclude_species );

			TS_ASSERT_EQUALS( cdr_options->include_only_germlines(), include_germlines );
			TS_ASSERT_EQUALS( cdr_options->exclude_germlines(), exclude_germlines );

			TS_ASSERT_EQUALS( cdr_options->include_only_pdbs(), include_pdbs );
			TS_ASSERT_EQUALS( cdr_options->exclude_pdbs(), exclude_pdbs );

		}

		for ( core::Size i = 1; i <= 5; ++i ) {
			CDRSetOptionsOP cdr_options = options[ i ];
			TS_ASSERT_EQUALS( cdr_options->load(), true );
		}

		TS_ASSERT_EQUALS( options[ l3 ]->load(), false);
		TS_ASSERT_EQUALS( options[ l1 ]->exclude_clusters(), exclude_l1_clusters );
		TS_ASSERT_EQUALS( options[ l1 ]->include_only_clusters(), include_l1_clusters );
	}

	void test_cdr_loading(){
		CDRSetOptionsParserOP parser( new CDRSetOptionsParser() );

		utility::vector1< CDRSetOptionsOP > options = parser->parse_options("protocols/antibody/database/cdr_set_instructions_load_test.txt");
		utility::vector1< CDRSetOptionsOP>  bogus_options = parser->parse_options("protocols/antibody/database/cdr_set_instructions.txt");

		//Turn off loading of everything but L1 for speed:
		options[ h1 ]->load(false);
		options[ h2 ]->load(false);
		options[ h3 ]->load(false);
		options[ l2 ]->load(false);
		options[ l3 ]->load(false);

		bogus_options[ h1 ]->load(false);
		bogus_options[ h2 ]->load(false);
		bogus_options[ h3 ]->load(false);
		bogus_options[ l2 ]->load(false);
		bogus_options[ l3 ]->load(false);

		AntibodyDatabaseManagerOP manager( new AntibodyDatabaseManager(ab_info, true /* Force use of north paper ab db */) );
		TS_ASSERT_THROWS_NOTHING(manager->load_cdr_poses(options, pose, false));

		TR << "------------ Checking to make sure an error is thrown when loading bogus options from database ------------- " << std::endl;
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(manager->load_cdr_poses(bogus_options, pose, false));
		TR << "------------ The preceeding error message was expected. -------------" << std::endl;


	}

};
