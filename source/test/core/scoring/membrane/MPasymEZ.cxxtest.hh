// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/membrane/MPasymEZ.cxxtest.hh
/// @brief  Test the scoring of the asymetric EZ potential
/// @author Meghan W Franklin

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/energy_methods/FaMPAsymEzCBEnergy.hh>
#include <core/energy_methods/FaMPAsymEzCGEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

static basic::Tracer TR("core.scoring.membrane.MPasymEZ.cxxtest");

using namespace core;

class MPasymEZTests : public CxxTest::TestSuite {

public:
	MPasymEZTests() {};

	void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_residue_energy()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::PoseOP pose_( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose_, "protocols/membrane/1C3W_TR_A.pdb", core::import_pose::PDB_file );

		//Add membrane information from spanfile:
		protocols::membrane::AddMembraneMoverOP add_memb( new protocols::membrane::AddMembraneMover( "protocols/membrane/1C3W_A.span" ) );
		add_memb->apply( *pose_ );


		core::Real const TOL(1e-3);
		//a failure of any pair of these sometimes indicates that the assigned Z-bin is one position off

		core::energy_methods::FaMPAsymEzCBEnergy asymCB;
		core::energy_methods::FaMPAsymEzCGEnergy asymCG;
		//check top of membrane (LEU9, CB_z = -13.559, CG_z = -14.139)
		EnergyMap emap;
		asymCB.residue_energy( pose_->residue( 9 ), *pose_, emap );
		asymCG.residue_energy( pose_->residue( 9 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCB ], -0.6398, TOL);
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCG ], -0.5787, TOL);

		//check top aromatic residue (TRP76, CB_z = -14.259, CG_z = -14.682)
		emap[ FaMPAsymEzCB ] = 0;
		emap[ FaMPAsymEzCG ] = 0;
		asymCB.residue_energy( pose_->residue( 76 ), *pose_, emap );
		asymCG.residue_energy( pose_->residue( 76 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCB ], -0.7776, TOL);
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCG ], -0.6557, TOL);

		//check middle membrane (ILE48, CB_z = -1.534, CG1_z = -0.613, CG2_z = -2.586)
		emap[ FaMPAsymEzCB ] = 0;
		emap[ FaMPAsymEzCG ] = 0;
		asymCB.residue_energy( pose_->residue( 48 ), *pose_, emap );
		asymCG.residue_energy( pose_->residue( 48 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCB ], -0.7200, TOL);
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCG ], -0.7100, TOL);

		//check Gly (GLY19, CA_z = 1.826)
		emap[ FaMPAsymEzCB ] = 0;
		emap[ FaMPAsymEzCG ] = 0;
		asymCB.residue_energy( pose_->residue( 19 ), *pose_, emap );
		asymCG.residue_energy( pose_->residue( 19 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCB ], -0.2000, TOL);
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCG ], -0.1900, TOL);

		//check Ala (ALA14, CB_z = -7.105)
		emap[ FaMPAsymEzCB ] = 0;
		emap[ FaMPAsymEzCG ] = 0;
		asymCB.residue_energy( pose_->residue( 14 ), *pose_, emap );
		asymCG.residue_energy( pose_->residue( 14 ), *pose_, emap );
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCB ], -0.3031, TOL);
		TS_ASSERT_DELTA( emap[ FaMPAsymEzCG ], -0.3160, TOL);

	}

};
