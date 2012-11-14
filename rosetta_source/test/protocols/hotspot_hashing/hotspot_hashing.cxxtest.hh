// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hotspot_hashing/hotspot_hashing.cxxtest.hh
/// @brief Test suite for protocols/hotspot_hashing/*
/// @author Alex Ford (fordas@uw.edu)
//
// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/MinimizingPatternSearch.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

static basic::Tracer TR("protocols.hotspot_hashing.hotspot_hashing.cxxtest");

namespace 
{
	/*
	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using core::conformation::ResidueOP;
	using core::conformation::ResidueFactory;
	using core::chemical::ChemicalManager;
	using core::chemical::ResidueType;
	*/

  class HotspotHashingTests : public CxxTest::TestSuite {
    public:
      void setUp() 
      {
        core_init();
      }

      void test_placeresidueattransform() 
      {
				using namespace protocols::hotspot_hashing;
        // Initialize target pose
        core::pose::Pose targetPose = create_test_in_pdb_pose();

        // Initialize residue representation
        core::conformation::ResidueOP residue;
				core::chemical::ResidueTypeSetCAP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
        core::chemical::ResidueType const & restype( residue_set->name_map( "ALA" ) );
        residue = core::conformation::ResidueFactory::create_residue( restype );

				Vector xunit = Vector(1, 0, 0);
				Vector yunit = Vector(0, 1, 0);
				Vector zunit = Vector(0, 0, 1);

				// Null transform
				{
					core::pose::Pose testPose(targetPose);
					TransformPair transform;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						TS_ASSERT_DELTA(placedResidue.xyz(placedResidue.atom_index("CA")), Vector(0, 0, 0), 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						Vector cb_vector(placedResidue.xyz(placedResidue.atom_index("CB")));
						TS_ASSERT_DELTA(cb_vector.normalize(), Vector(1, 0, 0), 1e-3);

						TS_ASSERT_DELTA(placedResidue.xyz(placedResidue.atom_index("C")).dot(Vector(0, 0, 1)), 0, 1e-3);
				}

				//Translate alone
				{
					Vector translation(10, 5, 2.5);

					core::pose::Pose testPose(targetPose);
					TransformPair transform;
					transform.translation = translation;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						TS_ASSERT_DELTA(placedResidue.xyz(placedResidue.atom_index("CA")), translation, 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						TS_ASSERT_DELTA((placedResidue.xyz("CB") - translation).normalize(), Vector(1, 0, 0), 1e-3);

						TS_ASSERT_DELTA((placedResidue.xyz("C") - translation).dot(Vector(0, 0, 1)), 0, 1e-3);
				}

				//Rotate in xy plane, centroid moved to y axis, CA in xy plane
				{
					Matrix rotation = rotation_matrix(xunit.cross(yunit), angle_of(xunit, yunit));

					core::pose::Pose testPose(targetPose);
					TransformPair transform;
					transform.rotation = rotation;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						TS_ASSERT_DELTA(placedResidue.xyz(placedResidue.atom_index("CA")), Vector(0, 0, 0), 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						TS_ASSERT_DELTA((placedResidue.xyz("CB") + Vector(0, 0, 0)).normalize(), Vector(0, 1, 0), 1e-3);

						TS_ASSERT_DELTA((placedResidue.xyz("C")).dot(Vector(0, 0, 1)), 0, 1e-3);
				}

				//Rotate in yz plane, centroid unmoved, CA in xz plane 
				{
					Matrix rotation = numeric::rotation_matrix(yunit.cross(zunit), angle_of(yunit, zunit));

					core::pose::Pose testPose(targetPose);
					TransformPair transform;
					transform.rotation = rotation;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						Vector ca_xyz = placedResidue.xyz("CA");
						TS_ASSERT_DELTA(ca_xyz, Vector(0, 0, 0), 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						Vector cb_xyz = placedResidue.xyz("CB");
						TS_ASSERT_DELTA(cb_xyz.normalize(), Vector(1, 0, 0), 1e-3);

						Vector c_xyz = placedResidue.xyz("C");
						TS_ASSERT_DELTA(c_xyz.dot(Vector(0, 1, 0)), 0, 1e-3);
				}

				//Rotate in xz plane, centroid moved to z axis, CA in yz plane 
				{
					Matrix rotation = numeric::rotation_matrix(xunit.cross(zunit), angle_of(xunit, zunit));

					core::pose::Pose testPose(targetPose);
					TransformPair transform;
					transform.rotation = rotation;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						Vector ca_xyz = placedResidue.xyz("CA");
						TS_ASSERT_DELTA(ca_xyz, Vector(0, 0, 0), 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						Vector cb_xyz = placedResidue.xyz("CB");
						TS_ASSERT_DELTA(cb_xyz.normalize(), Vector(0, 0, 1), 1e-3);

						Vector c_xyz = placedResidue.xyz("C");
						TS_ASSERT_DELTA(c_xyz.dot(Vector(1, 0, 0)), 0, 1e-3);
				}

				//Translate and rotate
				{
					Vector translation(10, 5, 2.5);
					Matrix rotation = numeric::rotation_matrix(xunit.cross(zunit), angle_of(xunit, zunit));

					core::pose::Pose testPose(targetPose);
					TransformPair transform;
					transform.translation = translation;
					transform.rotation = rotation;

					core::Size jumpindex;
					core::Size residueindex;

					MinimizingPatternSearch::placeResidueAtTransform(
							testPose,
							residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						TS_ASSERT_DELTA(placedResidue.xyz(placedResidue.atom_index("CA")), translation, 1e-3);

						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						TS_ASSERT_DELTA((placedResidue.xyz("CB") - translation).normalize(), Vector(0, 0, 1), 1e-3);

						TS_ASSERT_DELTA((placedResidue.xyz("C") - translation).dot(Vector(1, 0, 0)), 0, 1e-3);
				}
			}
  };
}
