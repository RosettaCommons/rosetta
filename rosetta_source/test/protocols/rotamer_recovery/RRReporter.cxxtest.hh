// -*- mode:c++;tab-width:2;indent-tabs-mode:nil;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RRReporter.cxxtest.hh
/// @brief  Reporter Classes for rotamer recovery
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/rotamer_recovery/RRReporter.hh>

// Project Headers
#include <test/core/init_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// C++ Headers
#include <iostream>

//Auto Headers
#include <utility/vector1.hh>



static basic::Tracer TR("protocols.rotamer_recovery.RRReporter.cxxtest");

class RRReporterTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
		core_init();
	}

	void test_RRReporterSimple_main() {
		do_test_RRReporterSimple_easy();

	}

	void
	do_test_RRReporterSimple_easy(){

		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using core::conformation::Residue;
		using protocols::rotamer_recovery::RRReporterSimple;

		RRReporterSimple rs;

		Pose pose ( fullatom_pose_from_string( pdb_string_1ten() ) );
		Residue residue ( pose.residue(1) );
		Real score(0);
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, true );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, true );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, false );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, true );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, false );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, false );
		rs.report_rotamer_recovery( pose, pose, residue, residue, score, true );

		rs.show( TR );

		TS_ASSERT_DELTA( rs.recovery_rate(), Real(4)/Real(7), .001 );
	}


};
