// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file test/core/scoring/cyclic_geometry_headers.cxxtest.hh
/// @brief Functions commonly used by cyclic peptide pose scoring unit tests.
/// @detials Cyclic permutations should score identically.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_test_core_scoring_cyclic_geometry_headers_HH
#define INCLUDED_test_core_scoring_cyclic_geometry_headers_HH

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>

using namespace std;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::AA;


static THREAD_LOCAL basic::Tracer TR_u("core.scoring.CyclicGeometryTests_utilfxns.cxxtest");

class CyclicGeometryTestHelper {

public:


/**********************************
FUNCTION DEFINITIONS
***********************************/

	/// @brief Test that the same score is returned for all cyclic permutations.
	///
	void cyclic_pose_test(
		core::scoring::ScoreFunctionOP sfxn,
		utility::vector1< core::pose::PoseOP > &poses,
		utility::vector1< core::pose::PoseOP > &mirror_poses
	) {
		//Score all of the poses and confirm that they're all equal to the first
		for ( core::Size i=1; i<=9; ++i ) {
			(*sfxn)(*poses[i]);
			if ( i>1 ) {
				TR_u << "\tTesting scoring with circular permutation of " << i - 1 << " residues." << std::endl;
				TS_ASSERT_DELTA(poses[1]->energies().total_energy(), poses[i]->energies().total_energy(), std::abs( std::max(poses[1]->energies().total_energy(), poses[i]->energies().total_energy())/10000.0 ) );
			}
			//Check mirrored geometry, too:
			TR_u << "\tTesting scoring with circular permutation of " << i - 1 << " residues and mirroring." << std::endl;
			(*sfxn)(*mirror_poses[i]);
			TS_ASSERT_DELTA(poses[1]->energies().total_energy(), mirror_poses[i]->energies().total_energy(), std::abs( std::max(poses[1]->energies().total_energy(), mirror_poses[i]->energies().total_energy())/10000.0 ) );
			TR_u.flush();
		}
	}

}; //class CyclicGeometryTestHelper


#endif //INCLUDED_test_core_scoring_cyclic_geometry_headers_HH
