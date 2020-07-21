// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/chemical/NCAACentroidResidueTypeTests.cxxtest.hh
/// @brief  Unit tests confirming that centroid resiues for NCAA types are generated reasonably.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers

// Protocols Headers (for convenience)
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("NCAACentroidResidueTypeTests");


class NCAACentroidResidueTypeTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-symmetric_gly_tables true" );
	}

	void tearDown(){

	}

	void test_centroid_aib_orn_dab_dap(){
		TR << "Starting NCAACentroidResidueTypeTests::test_centroid_aib_orn_dab_dap()." << std::endl;
		TR << "This test confirms that we can build poses with AIB, ORN, DAP, and DAB, that we can convert these types to centroid representations, that we can mirror the centroid representations to get the D-equivalents, and that we can convert back to full atom representations." << std::endl;

		using namespace protocols::cyclic_peptide;
		using namespace protocols::simple_moves;

		PeptideStubMover stubmover;
		stubmover.set_reset_mode(true);
		stubmover.add_residue( "Append", "GLY", 1, true, "", 0, 0, nullptr, "" );
		stubmover.add_residue( "Append", "ORN", 1, false, "", 0, 0, nullptr, "" );
		stubmover.add_residue( "Append", "DAB", 1, false, "", 0, 0, nullptr, "" );
		stubmover.add_residue( "Append", "DPP", 1, false, "", 0, 0, nullptr, "" );
		stubmover.add_residue( "Append", "AIB", 1, false, "", 0, 0, nullptr, "" );
		stubmover.add_residue( "Append", "GLY", 1, false, "", 0, 0, nullptr, "" );

		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		stubmover.apply(*pose);

		TS_ASSERT_EQUALS( pose->total_residue(), 6 );
		for ( core::Size i(1), imax(pose->total_residue()); i<=imax; ++i ) {
			pose->set_phi( i, -61 );
			pose->set_psi( i, -41 );
			pose->set_omega( i, 180 );
		}
		pose->update_residue_neighbors();
		//pose->dump_pdb( "VTEMP1.pdb" );

		SwitchResidueTypeSetMover switcher;
		switcher.type_set_tag( "centroid" );
		TS_ASSERT_THROWS_NOTHING( switcher.apply(*pose) );
		//pose->dump_pdb( "VTEMP2.pdb" );
		TS_ASSERT_EQUALS( pose->total_residue(), 6 );
		TS_ASSERT_EQUALS( pose->residue_type(1).base_name(), "GLY" );
		TS_ASSERT_EQUALS( pose->residue_type(2).base_name(), "ORN" );
		TS_ASSERT_EQUALS( pose->residue_type(3).base_name(), "DAB" );
		TS_ASSERT_EQUALS( pose->residue_type(4).base_name(), "DPP" );
		TS_ASSERT_EQUALS( pose->residue_type(5).base_name(), "AIB" );
		TS_ASSERT_EQUALS( pose->residue_type(6).base_name(), "GLY" );

		FlipChiralityMover flipper;
		TS_ASSERT_THROWS_NOTHING( flipper.apply( *pose ) );
		//pose->dump_pdb( "VTEMP3.pdb" );
		TS_ASSERT_EQUALS( pose->total_residue(), 6 );
		TS_ASSERT_EQUALS( pose->residue_type(1).base_name(), "GLY" );
		TS_ASSERT_EQUALS( pose->residue_type(2).base_name(), "DORN" );
		TS_ASSERT_EQUALS( pose->residue_type(3).base_name(), "DDAB" );
		TS_ASSERT_EQUALS( pose->residue_type(4).base_name(), "DDPP" );
		TS_ASSERT_EQUALS( pose->residue_type(5).base_name(), "AIB" );
		TS_ASSERT_EQUALS( pose->residue_type(6).base_name(), "GLY" );

		SwitchResidueTypeSetMover switcher2;
		TS_ASSERT_THROWS_NOTHING( switcher2.type_set_tag("fa_standard") );
		//pose->dump_pdb( "VTEMP4.pdb" );
		TS_ASSERT_EQUALS( pose->total_residue(), 6 );
		TS_ASSERT_EQUALS( pose->residue_type(1).base_name(), "GLY" );
		TS_ASSERT_EQUALS( pose->residue_type(2).base_name(), "DORN" );
		TS_ASSERT_EQUALS( pose->residue_type(3).base_name(), "DDAB" );
		TS_ASSERT_EQUALS( pose->residue_type(4).base_name(), "DDPP" );
		TS_ASSERT_EQUALS( pose->residue_type(5).base_name(), "AIB" );
		TS_ASSERT_EQUALS( pose->residue_type(6).base_name(), "GLY" );

		TR << "Finished NCAACentroidResidueTypeTests::test_centroid_aib_orn_dab_dap()." << std::endl;
	}

};
