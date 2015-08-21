// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/fiber_diffraction/FiberDiffractionEnergy.cxxtest.hh
/// @brief test suite for core::scoring::fiber_diffraction::FiberDiffractionEnergy
/// @author Wojciech Potrzebowski and Ingemar Andre

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <core/scoring/fiber_diffraction/FiberDiffractionEnergy.hh>

// Project headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <sstream>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/EXCN_Base.hh>


static basic::Tracer TR("test.scoring.fiber_diffraction.FiberDiffractionEnergy");

using namespace core;
using namespace scoring;
using namespace scoring::symmetry;
using namespace pose;


class FiberDiffractionEnergyTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}


	void tearDown() {}

	// @brief test default options and default locator
	void test_FiberDiffractionEnergy() {

		try{
			do_test_FiberDiffractionEnergy();
			do_test_FiberDiffractionEnergyDens();
		} catch( utility::excn::EXCN_Base& excn ) {
			TR << "ERROR: Exception caught by FiberDiffractionEnergy:"
				<< excn << std::endl;
		}
	}
	void do_test_FiberDiffractionEnergy() {

		core_init_with_additional_options( "-fiber_diffraction::layer_lines core/scoring/fiber_diffraction/data/1ifp.mini.dat -fiber_diffraction:a 27 -fiber_diffraction:b 5 -fiber_diffraction:resolution_cutoff_low 0.0833333 -fiber_diffraction:resolution_cutoff_high 0.333333 -symmetry:symmetry_definition core/scoring/fiber_diffraction/data/helix_denovo.sdef" );


		std::string pdb_file_name = "core/scoring/fiber_diffraction/data/1IFP.mini.pdb";
		Pose pose;
		core::import_pose::pose_from_pdb(pose, pdb_file_name);
		core::pose::symmetry::make_symmetric_pose( pose );

		core::Real const TOL(1e-5);

		SymmetricScoreFunction sfxn;
		sfxn.set_weight(fiberdiffraction, 1 );
		Energy score = sfxn( pose );

		TS_ASSERT( std::fabs( score - 0.564044 ) < TOL );

	}

	void do_test_FiberDiffractionEnergyDens() {

		core_init_with_additional_options( "-fiber_diffraction::layer_lines core/scoring/fiber_diffraction/data/1ifp.mini.dat -in:file:centroid_input -fiber_diffraction:a 27 -fiber_diffraction:b 5 -fiber_diffraction:resolution_cutoff_low 0.0833333 -fiber_diffraction:resolution_cutoff_high 0.333333 -symmetry:symmetry_definition core/scoring/fiber_diffraction/data/helix_denovo.sdef " );

		std::string pdb_file_name = "core/scoring/fiber_diffraction/data/1IFP.mini.pdb";
		Pose pose;
		core::import_pose::pose_from_pdb(pose, pdb_file_name);
		core::pose::symmetry::make_symmetric_pose( pose );

		core::Real const TOL(1e-5);

		SymmetricScoreFunction sfxn;
		sfxn.set_weight(fiberdiffractiondens, 1 );

		Energy score = sfxn( pose );

		TS_ASSERT( std::fabs( score -  0.99333 ) < TOL );
	}
};
