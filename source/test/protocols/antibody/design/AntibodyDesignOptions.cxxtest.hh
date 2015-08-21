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

#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#define BFE BOOST_FOREACH
using namespace protocols::antibody;
using namespace protocols::antibody::design;
using utility::vector1;

static thread_local basic::Tracer TR("protocols.antibody.AntibodyConstraintTests");
class AntibodyDesignOptions: public CxxTest::TestSuite {
	core::pose::Pose pose;
	AntibodyInfoOP ab_info;


public:

	void setUp(){
		core_init();
		core::import_pose::pose_from_pdb(pose, "protocols/antibody/aho_with_antigen.pdb"); //AHO renumbered pose
		ab_info = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );

	}

	void test_options_seq_design() {
		CDRSeqDesignOptionsParserOP parser( new CDRSeqDesignOptionsParser() );
		utility::vector1< CDRSeqDesignOptionsOP > options = parser->parse_options("protocols/antibody/design/seq_des_instructions.txt");
		for ( core::Size i = 1; i <= 5; ++i ) {
			CDRSeqDesignOptionsOP cdr_option = options[ i ];
			TS_ASSERT_EQUALS( cdr_option->design(), true );
		}

		for ( core::Size i = 1; i <= 4; ++i ) {
			TS_ASSERT_EQUALS( options[ i ]->design_strategy(), seq_design_profiles );
		}

		TS_ASSERT_EQUALS( options[ l3 ]->design(), false );

		TS_ASSERT_EQUALS( options[ l2 ]->design_strategy(), seq_design_conservative );
		TS_ASSERT_EQUALS( options[ l3 ]->design_strategy(), seq_design_basic );

	}

	void test_options_graft_design() {
		CDRGraftDesignOptionsParserOP parser( new CDRGraftDesignOptionsParser() );
		utility::vector1< CDRGraftDesignOptionsOP > options = parser->parse_options("protocols/antibody/design/graft_des_instructions.txt");
		for ( core::Size i = 1; i <=6; ++i ) {
			CDRGraftDesignOptionsOP cdr_option = options[ i ];
			TS_ASSERT_EQUALS( cdr_option->mintype(), relax );
		}

		for ( core::Size i = 1; i <=5; ++i ) {
			CDRGraftDesignOptionsOP cdr_option = options[ i ];
			TS_ASSERT_EQUALS( cdr_option->design(), true );
		}

		//Test Bogus L1 settings:
		CDRGraftDesignOptionsOP l3_option = options[ l3 ];
		TS_ASSERT_EQUALS( l3_option->design(), false );
		TS_ASSERT_DELTA( l3_option->weight(), 2.0, .0001 );

	}

};
