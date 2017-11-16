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


static basic::Tracer TR_u("core.scoring.CyclicGeometryTests_utilfxns.cxxtest");

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
		utility::vector1< core::pose::PoseOP > &mirror_poses,
		bool const do_cutpoint_test
	) {
		core::Size nposes(poses.size());
		runtime_assert(nposes == mirror_poses.size());

		//Score all of the poses and confirm that they're all equal to the first
		for ( core::Size i=1; i<=nposes; ++i ) {
			(*sfxn)(*poses[i]);
			core::Real const energy1( poses[1]->energies().total_energy() );
			if ( i>1 ) {
				TR_u << "\tTesting scoring with circular permutation of " << i - 1 << " residues:\t";
				core::Real const energy2( poses[i]->energies().total_energy() );
				TR_u << "E1=" << energy1 << "\tE" << i-1 << "=" << energy2 << std::endl;
				TS_ASSERT_DELTA(energy1, energy2, std::abs( std::max(energy1, energy2)/10000.0 ) );
			}
			//Check mirrored geometry, too:
			TR_u << "\tTesting scoring with circular permutation of " << i - 1 << " residues and mirroring:\t";
			(*sfxn)(*mirror_poses[i]);
			core::Real const energy3( mirror_poses[i]->energies().total_energy() );
			TR_u << "E1=" << energy1 << "\tmirror_E" << i-1 << "=" << energy3 << std::endl;
			TS_ASSERT_DELTA(energy1, energy3, std::abs( std::max(energy1, energy3)/10000.0 ) );
			TR_u.flush();
		}

		if ( do_cutpoint_test ) {
			//Add the cutpoint variants to all the poses, score, and confirm that they're equal again.
			core::pose::PoseOP first_cutpoint_pose;
			for ( core::Size i(1); i<=nposes; ++i ) {
				core::pose::PoseOP cutpoint_i( poses[i]->clone() );
				core::pose::PoseOP mirr_cutpoint_i( mirror_poses[i]->clone() );
				core::pose::correctly_add_cutpoint_variants( *cutpoint_i, nposes, false, 1 );
				core::pose::correctly_add_cutpoint_variants( *mirr_cutpoint_i,nposes, false, 1 );
				cutpoint_i->set_phi(1, poses[i]->phi(1));
				cutpoint_i->set_psi(nposes, poses[i]->psi(nposes));
				cutpoint_i->set_omega(nposes, poses[i]->omega(nposes) );
				mirr_cutpoint_i->set_phi(1, mirror_poses[i]->phi(1));
				mirr_cutpoint_i->set_psi(nposes, mirror_poses[i]->psi(nposes));
				mirr_cutpoint_i->set_omega(nposes, mirror_poses[i]->omega(nposes) );
				(*sfxn)(*cutpoint_i);
				(*sfxn)(*mirr_cutpoint_i);
				if ( i==1 ) {
					first_cutpoint_pose = cutpoint_i;
				} else {
					TR_u << "\tTesting scoring with cutpoints and circular permutation of " << i-1 << " residues:\t";
					core::Real const energy2( cutpoint_i->energies().total_energy() );
					core::Real const energy1( first_cutpoint_pose->energies().total_energy() );
					TR_u << "E1=" << energy1 << "\tE" << i-1 << "=" << energy2 << std::endl;
					TS_ASSERT_DELTA( energy2, energy1, std::abs( std::max( energy2, energy1 )/10000.0 ) );
				}
				TR_u << "\tTesting scoring with cutpoints and circular permutation of " << i-1 << " residues, with mirroring:\t";
				core::Real const energy1( first_cutpoint_pose->energies().total_energy() );
				core::Real const energy3( mirr_cutpoint_i->energies().total_energy() );
				TR_u << "E1=" << energy1 << "\tE" << i-1 << "=" << energy3 << std::endl;
				TS_ASSERT_DELTA( energy3, energy1, std::abs( std::max( energy3, energy1 )/10000.0 ) );

				//The following is for debugging only, and can be deleted:
				/*
				char fname[256]; //DELETE ME
				sprintf(fname, "mirr_%04lu.pdb", i - 1); //DELETE ME
				mirr_cutpoint_i->dump_scored_pdb( std::string( fname ), (*sfxn) ); //DELETE ME
				sprintf(fname, "ref_%04lu.pdb", i - 1); //DELETE ME
				mirror_poses[i]->dump_scored_pdb( std::string(fname), (*sfxn) ); //DELETE ME
				core::Size const nres(mirr_cutpoint_i->total_residue());
				TR_u << "PHI1=" << mirr_cutpoint_i->phi(1) << "\tPHI1ref=" << mirror_poses[i]->phi(1) << std::endl; //DELETE ME
				TR_u << "PSI1=" << mirr_cutpoint_i->psi(1) << "\tPSI1ref=" << mirror_poses[i]->psi(1) << std::endl; //DELETE ME
				TR_u << "OMG1=" << mirr_cutpoint_i->omega(1) << "\tOMG1ref=" << mirror_poses[i]->omega(1) << std::endl; //DELETE ME
				TR_u << "PHIN=" << mirr_cutpoint_i->phi(nres) << "\tPHINref=" << mirror_poses[i]->phi(nres) << std::endl; //DELETE ME
				TR_u << "PSIN=" << mirr_cutpoint_i->psi(nres) << "\tPSINref=" << mirror_poses[i]->psi(nres) << std::endl; //DELETE ME
				TR_u << "OMGN=" << mirr_cutpoint_i->omega(nres) << "\tOMGNref=" << mirror_poses[i]->omega(nres) << std::endl; //DELETE ME
				*/
			}
		}
	}

	void
	cyclic_pose_test(
		core::scoring::ScoreFunctionOP sfxn,
		utility::vector1< core::pose::PoseOP > &poses,
		utility::vector1< core::pose::PoseOP > &mirror_poses
	) {
		cyclic_pose_test( sfxn, poses, mirror_poses, true );
	}

}; //class CyclicGeometryTestHelper


#endif //INCLUDED_test_core_scoring_cyclic_geometry_headers_HH
