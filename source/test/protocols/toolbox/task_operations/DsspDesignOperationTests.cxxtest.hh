// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @author Brian Koepnick (koepnick@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/toolbox/task_operations/DsspDesignOperation.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <string>
#include <map>

using namespace protocols::toolbox::task_operations;

class DsspDesignOperationTests : public CxxTest::TestSuite {

public:

	typedef std::map< std::string, std::string > SecStructResidues;

public:

	void setUp(){
		protocols_init();
	}

	void test_DsspDesignOperation() {
		// check output of DsspDesignOperation Tracer
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/toolbox/task_operations/dssp_in.pdb" , core::import_pose::PDB_file);

		basic::otstreamOP UT( new test::UTracer( "protocols/toolbox/task_operations/DsspDesignOperation.u" ) );
		basic::Tracer::set_ios_hook( UT, "protocols.toolbox.TaskOperations.DsspDesignOperation" );

		TaskFactory Dssp_factory;
		Dssp_factory.push_back( core::pack::task::operation::TaskOperationCOP( new DsspDesignOperation ) );
		Dssp_factory.create_task_and_apply_taskoperations( pose );

		basic::Tracer::set_ios_hook( 0, "" );
	}

	void test_copy_constructor() {

		protocols::jd2::parser::BluePrintOP blueprint;
		blueprint = protocols::jd2::parser::BluePrintOP( new protocols::jd2::parser::BluePrint(
			"protocols/toolbox/task_operations/dssp_in.blueprint" ) );

		DsspDesignOperationOP dssp_design( new DsspDesignOperation );
		dssp_design->set_blueprint( blueprint );
		dssp_design->set_restrictions_aa( "Nterm", "Y" );

		DsspDesignOperationOP clone( new DsspDesignOperation( *dssp_design ) );

		TS_ASSERT_EQUALS( dssp_design->blueprint_, clone->blueprint_ );
		TS_ASSERT_EQUALS( clone->sse_residues_[ "Nterm" ], "Y" )
			}

			void test_get_restrictions() {
			DsspDesignOperationOP dssp_design( new DsspDesignOperation );
		utility::vector1< bool > helixcap_res = dssp_design->get_restrictions( "HelixCapping" );

		TS_ASSERT( helixcap_res.size() == 20 );

		// default HelixCapping residues are DNST
		for ( core::Size aa=1; aa<=20; ++aa ) {
			if ( aa == 3 || aa == 12 || aa == 16 || aa == 17 ) {
				TS_ASSERT( helixcap_res[ aa ] );
			} else {
				TS_ASSERT( !helixcap_res[ aa ] );
			}
		}
	}

	void test_set_restrictions_aa() {
		// set allowed residue for one SSE
		DsspDesignOperationOP strand_gly( new DsspDesignOperation );
		strand_gly->set_restrictions_aa( "Strand", "G" );

		TS_ASSERT_EQUALS( strand_gly->sse_residues_[ "Strand" ], "G" );
		TS_ASSERT_DIFFERS( strand_gly->sse_residues_[ "Helix" ], "G" );

		// set allowed residues for all SSEs
		DsspDesignOperationOP all_pg( new DsspDesignOperation );
		all_pg->set_restrictions_aa( "all", "PG" );

		for ( SecStructResidues::iterator it = all_pg->sse_residues_.begin(); it != all_pg->sse_residues_.end(); it++ ) {
			TS_ASSERT_EQUALS( it->second, "PG" );
		}
	}

	void test_set_restrictions_append() {
		// append allowed residue for one SSE
		DsspDesignOperationOP helix_pro( new DsspDesignOperation );
		helix_pro->set_restrictions_append( "Helix", "P" );

		TS_ASSERT_DIFFERS( helix_pro->sse_residues_[ "Helix" ].find( "P" ), std::string::npos );
		TS_ASSERT_EQUALS( helix_pro->sse_residues_[ "Strand" ].find( "P" ), std::string::npos );

		// append allowed residues for all SSEs
		DsspDesignOperationOP all_pg( new DsspDesignOperation );
		all_pg->set_restrictions_append( "all", "PG" );

		for ( SecStructResidues::iterator it = all_pg->sse_residues_.begin(); it != all_pg->sse_residues_.end(); it++ ) {
			TS_ASSERT_DIFFERS( it->second.find( "P" ), std::string::npos );
			TS_ASSERT_DIFFERS( it->second.find( "G" ), std::string::npos );
		}
	}

	void test_set_restrictions_exclude() {
		// exclude allowed residue for one SSE
		DsspDesignOperationOP loop_nopro( new DsspDesignOperation );
		loop_nopro->set_restrictions_exclude( "Loop", "P" );

		TS_ASSERT_EQUALS( loop_nopro->sse_residues_[ "Loop" ].find( "P" ), std::string::npos )
			TS_ASSERT_DIFFERS( loop_nopro->sse_residues_[ "HelixStart" ].find( "P" ), std::string::npos )

			// append allowed residues for all SSEs
			DsspDesignOperationOP all_st( new DsspDesignOperation );
		all_st->set_restrictions_exclude( "all", "ST" );

		for ( SecStructResidues::iterator it = all_st->sse_residues_.begin(); it != all_st->sse_residues_.end(); it++ ) {
			TS_ASSERT_EQUALS( it->second.find( "S" ), std::string::npos );
			TS_ASSERT_EQUALS( it->second.find( "T" ), std::string::npos );
		}
	}

	// it'd be nice to test if the XML Schema itself is valid.
	void dont_test_dssp_design_operation_xsd() {
		using namespace utility::tag;
		XMLSchemaDefinition xsd;
		DsspDesignOperation::provide_xml_schema( xsd );
		std::cout << "DsspDesignOperation XSD:\n" << xsd.full_definition() << std::endl;
	}

};
