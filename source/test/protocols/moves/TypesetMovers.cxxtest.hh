// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/TypesetMovers.cxxtest.hh
/// @brief  test for ReturnSidechainMover and SwitchToResidueTypesetMover
/// @author Steven Lewis (smlewi@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Package headers
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

class TypesetMoversTests : public CxxTest::TestSuite {

	core::pose::Pose pose;

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		core::import_pose::pose_from_file( pose, "protocols/moves/test_in_short.pdb" , core::import_pose::PDB_file);
	}

	void tearDown() {
		pose.clear();
	}


	// --------------- Test Cases --------------- //

	/// @brief take a pose, convert it to centroid (and check the pose), run ReturnSidechainMover to return its sidechains (and check the pose again)
	void test_TypesetMovers() {

		core::pose::Pose const posecopy(pose);
		protocols::simple_moves::ReturnSidechainMover RSmover(pose);
		protocols::simple_moves::SwitchResidueTypeSetMover res_mover(core::chemical::CENTROID);

		res_mover.apply( pose ); //now in centroid mode
		//check that the pose was centroidized correctly
		test::UTracer UT_SRTSM("protocols/moves/centroid_typeset.test_in_short.pdb");
		pose.dump_pdb(UT_SRTSM);

		RSmover.apply( pose ); //now re-fullatom-ized, should retain proper sidechains
		//check that the pose was returned to fullatom correctly
		test::UTracer UT_RSM("protocols/moves/fullatom_typeset.test_in_short.pdb");
		pose.dump_pdb(UT_RSM);

		//test another way that it came out right: for all residues, for all chi, check equality
		core::Size nres = pose.total_residue();
		TS_ASSERT_EQUALS(nres, posecopy.total_residue());
		for ( core::Size res(1); res <= nres; ++res ) {
			utility::vector1<core::Real> const & copychi(posecopy.residue(res).chi());
			utility::vector1<core::Real> const & newchi(pose.residue(res).chi());
			TS_ASSERT_EQUALS(copychi.size(), newchi.size());
			for ( core::Size chi(1); chi <= copychi.size(); ++chi ) {
				TS_ASSERT_DELTA(copychi.at(chi), newchi.at(chi), 0.0001);
			}
		}
	}//end test_TypesetMovers

};//end class
