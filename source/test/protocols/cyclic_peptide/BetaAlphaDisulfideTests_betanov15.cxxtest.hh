// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/BetaAlphaDisulfideTests_betanov15.cxxtest.hh
/// @brief  Unit tests for working with peptides mixed alpha-cys/beta-cys disulfides.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#define EXPECTEDDIST 2.01
#define DISTTHRESH 0.05
#define EXPECTEDANGLE 1.81898214643 //104.22 degrees
#define ANGTHRESH 0.01745329251 // 1 degree


static basic::Tracer TR("protocols.cyclic_peptide.BetaAlphaDisulfide_betanov15_Tests");

class BetaAlphaDisulfide_betanov15_Tests : public CxxTest::TestSuite {

public:

	void setUp(){
		core_init_with_additional_options("-write_all_connect_info");
	}

	void tearDown(){
	}

	/// @brief Score a structure involving mixed beta-cys -- alpha-cys disulfides.
	///
	void test_mixed_disulf() {
		core::pose::PoseOP pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose, "protocols/cyclic_peptide/mixed_disulf_test.pdb" , core::import_pose::PDB_file);

		//Get rid of other beta-amino acids, which can't repack because I don't distribute rotamer libraries for all of them with Rosetta:
		protocols::simple_moves::MutateResidue mutres33( 33, "B3A" );
		protocols::simple_moves::MutateResidue mutres35( 35, "B3A" );
		protocols::simple_moves::MutateResidue mutres39( 39, "B3A" );
		mutres33.apply(*pose);
		mutres35.apply(*pose);
		mutres39.apply(*pose);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::dslf_fa13, 1.0 );

		protocols::relax::FastRelax frlx(sfxn, 3);
		frlx.apply( *pose );

		core::Real dist1, dist2;
		dist1 = pose->residue(32).xyz("SG").distance( pose->residue(3).xyz("SG") );
		dist2 = pose->residue(36).xyz("SG").distance( pose->residue(22).xyz("SG") );
		TR << "SG-SG distance, expected (A): " << EXPECTEDDIST << std::endl;
		TR << "SG-SG distance, beta32-dalpha3 (A): " << dist1 << std::endl;
		TR << "SG-SG distance, beta36-alpha22 (A): " << dist2 << std::endl;
		TS_ASSERT_DELTA( dist1, EXPECTEDDIST, DISTTHRESH);
		TS_ASSERT_DELTA( dist2, EXPECTEDDIST, DISTTHRESH);

		core::Real ang1, ang2, ang3, ang4;
		ang1 = numeric::angle_of( pose->residue(32).xyz("CB"), pose->residue(32).xyz("SG"), pose->residue(3).xyz("SG") );
		ang2 = numeric::angle_of( pose->residue(32).xyz("SG"), pose->residue(3).xyz("SG"), pose->residue(3).xyz("CB") );
		ang3 = numeric::angle_of( pose->residue(36).xyz("CB"), pose->residue(36).xyz("SG"), pose->residue(22).xyz("SG") );
		ang4 = numeric::angle_of( pose->residue(36).xyz("SG"), pose->residue(22).xyz("SG"), pose->residue(22).xyz("CB") );
		TR << "CB-SG-SG angle, expected (rad): " << EXPECTEDANGLE << std::endl;
		TR << "CB-SG-SG angle, beta32 (rad): " << ang1 << std::endl;
		TR << "CB-SG-SG angle, dalpha3 (rad): " << ang2 << std::endl;
		TR << "CB-SG-SG angle, beta36 (rad): " << ang3 << std::endl;
		TR << "CB-SG-SG angle, alpha22 (rad): " << ang4 << std::endl;
		TS_ASSERT_DELTA( ang1, EXPECTEDANGLE, ANGTHRESH);
		TS_ASSERT_DELTA( ang2, EXPECTEDANGLE, ANGTHRESH);
		TS_ASSERT_DELTA( ang3, EXPECTEDANGLE, ANGTHRESH);
		TS_ASSERT_DELTA( ang4, EXPECTEDANGLE, ANGTHRESH);
	}

private:

};



