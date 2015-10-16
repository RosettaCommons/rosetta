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
#include <protocols/antibody/design/NativeAntibodySeq.hh>
#include <protocols/antibody/design/util.hh>

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

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.AntibodyConstraintTests");
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
		TS_ASSERT_EQUALS( options[ l3 ]->fallback_strategy(), seq_design_none );

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
	
	void test_native_seq() {
		NativeAntibodySeq native_record = NativeAntibodySeq(pose, ab_info);
		TS_ASSERT_THROWS_NOTHING(native_record.set_sequence(pose) ); //Does nothing.  Just a test.
		TS_ASSERT_THROWS_NOTHING( native_record.set_to_pose( pose ) );
		TS_ASSERT_EQUALS( protocols::antibody::design::has_native_sequence( pose ), true);
		
		
		std::string pose_seq = pose.sequence();
		
		for (core::Size i = 1; i <= 6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			TS_ASSERT_THROWS_NOTHING(native_record.set_from_cdr(pose, cdr) ); //Does nothing. Just a test.
			TS_ASSERT_THROWS_NOTHING( set_native_cdr_sequence( ab_info, cdr, pose));
		}
		std::string record_seq = native_record.get_sequence(pose);
		TS_ASSERT_EQUALS( pose_seq, record_seq);
		
		std::string seq_from_pose_cache = protocols::antibody::design::get_native_sequence(pose);
		TS_ASSERT_EQUALS( pose_seq, seq_from_pose_cache);
		
		
		
	}

};
