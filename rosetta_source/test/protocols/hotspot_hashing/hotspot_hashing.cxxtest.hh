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
#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/RT.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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
	using core::kinematics::RT;

  class HotspotHashingTests : public CxxTest::TestSuite {

    public:
      void setUp() 
      {
        core_init();
      }

			void test_stubcentroid_transform()
			{
				using namespace protocols::hotspot_hashing;

        // Initialize residue representation
        core::conformation::ResidueOP residue;
				core::chemical::ResidueTypeSetCAP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
				
				Vector xunit = Vector(1, 0, 0);
				Vector yunit = Vector(0, 1, 0);
				Vector zunit = Vector(0, 0, 1);

				{
					core::chemical::ResidueType const & restype( residue_set->name_map( "ALA" ) );
					residue = core::conformation::ResidueFactory::create_residue( restype );

					Vector stub_centroid = PlaceMinimizeSearch::residueStubCentroid(*residue);
					RT centroid_transform = PlaceMinimizeSearch::residueStubCentroidTransform(*residue);


					Vector ca_location = residue->xyz("CA");
					// Using alanine residue, so centroid of SC heavyatoms will be CB position
					Vector cb_location = residue->xyz("CB");
					Vector c_location = residue->xyz("C");

					TS_ASSERT_DELTA(cb_location, stub_centroid, 1e-3);
					TS_ASSERT_DELTA(-ca_location, centroid_transform.get_translation(), 1e-3);

					Vector cb_vector = cb_location - ca_location;
					Vector cb_rotated = centroid_transform.get_rotation() * cb_vector;
					TS_ASSERT_DELTA(xunit, cb_rotated.normalize(), 1e-3);
				}
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
					RT transform;

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
							transform,
							jumpindex,
							residueindex);

						TS_ASSERT_EQUALS(jumpindex, testPose.num_jump());
						TS_ASSERT_EQUALS(residueindex, testPose.total_residue());

						TS_ASSERT_EQUALS(jumpindex, targetPose.num_jump() + 1);
						TS_ASSERT_EQUALS(residueindex, targetPose.total_residue() + 1);

						core::conformation::Residue const &placedResidue = testPose.residue(residueindex);

						TS_ASSERT_EQUALS(placedResidue.type().name(), "ALA");

						Vector ca_location = placedResidue.xyz("CA");
						// Using alanine residue, so centroid of SC heavyatoms will be CB position
						Vector cb_vector = placedResidue.xyz("CB");
						cb_vector = cb_vector.normalize();
						Vector c_location = placedResidue.xyz("C");

						TS_ASSERT_DELTA(ca_location, Vector(0, 0, 0), 1e-3);

						TS_ASSERT_DELTA(cb_vector, Vector(1, 0, 0), 1e-3);

						TS_ASSERT_DELTA(c_location.dot(Vector(0, 0, 1)), 0, 1e-3);
				}

				//Translate alone
				{
					Vector translation(10, 5, 2.5);

					core::pose::Pose testPose(targetPose);
					RT transform;
					transform.set_translation(translation);

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
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
					RT transform;
					transform.set_rotation(rotation);

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
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
					RT transform;
					transform.set_rotation(rotation);

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
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
					RT transform;
					transform.set_rotation(rotation);

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
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
					RT transform;
					transform.set_translation(translation);
					transform.set_rotation(rotation);

					core::Size jumpindex;
					core::Size residueindex;

					PlaceMinimizeSearch::placeResidueAtTransform(
							testPose,
							*residue,
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

      void test_lsmsearchpattern() 
      {
				using namespace protocols::hotspot_hashing;

				core::Real angle_sampling = 90;
				core::Real translocation_sampling = 1;
				core::Real max_radius = 2;
				core::Real distance_sampling = 1;
				core::Real max_distance = 2;


				VectorPair lsm_zero(
							Vector(0, 0, 0),
							Vector(1, 1, 1));
				LSMSearchPattern zero_pattern(
						lsm_zero,
						angle_sampling,
						translocation_sampling,
						max_radius,
						distance_sampling,
						max_distance);

				utility::vector1<RT> points_zero = zero_pattern.Searchpoints();
				TS_ASSERT_EQUALS(points_zero.size(), 4 /*angles*/ * 13 /*trans*/ * 3 /*dist*/);

				// Assert on location
				foreach(RT transform, points_zero)
				{
					TS_ASSERT_LESS_THAN_EQUALS(lsm_zero.position.distance(transform.get_translation()), std::sqrt( max_radius * max_radius + max_distance * max_distance));
				}

				//LSM taken from test structure
				VectorPair lsm_test(
						Vector(64.856, 29.883, 42.865),
						Vector(-0.56293, 0.20554, -0.80054));
				LSMSearchPattern test_pattern(
						lsm_test,
						90,
						1,
						2,
						1,
						2);
				utility::vector1<RT> points_test = test_pattern.Searchpoints();
				TS_ASSERT_EQUALS(points_test.size(), 4 /*angles*/ * 13 /*trans*/ * 3 /*dist*/);

				foreach(RT transform, points_test)
				{
					TS_ASSERT_LESS_THAN_EQUALS(lsm_test.position.distance(transform.get_translation()), std::sqrt( max_radius * max_radius + max_distance * max_distance));
				}
			}
			 
			void place_and_assert_transform(core::pose::Pose & targetPose, core::conformation::ResidueOP residue, protocols::hotspot_hashing::RT transform)
			{
				using namespace protocols::hotspot_hashing;
				core::Size jumpindex;
				core::Size residueindex;

				PlaceMinimizeSearch::placeResidueAtTransform(
						targetPose,
						*residue,
						transform,
						jumpindex,
						residueindex);

				core::conformation::Residue const &placedResidue = targetPose.residue(residueindex);

				Vector residue_CA_location = placedResidue.xyz(placedResidue.atom_index("CA"));
				TS_ASSERT_DELTA(residue_CA_location, transform.get_translation(), 1e-3);

				Vector stubcentroid = PlaceMinimizeSearch::residueStubCentroid(placedResidue);
				Vector CA_centroid_vector = (stubcentroid - residue_CA_location).normalize();
				Vector xunit_rotated = transform.get_rotation() * Vector(1, 0, 0);

				TS_ASSERT_DELTA(CA_centroid_vector, xunit_rotated, 1e-3); 		
			}

      void test_lsmsearchpattern_placement() 
			{
				using namespace protocols::hotspot_hashing;
        // Initialize target pose
        core::pose::Pose targetPose;
				core::import_pose::pose_from_pdb( targetPose, "protocols/hotspot_hashing/3ve0_IJ.pdb" );

        // Initialize residue representation
        core::conformation::ResidueOP residue;
				core::chemical::ResidueTypeSetCAP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
        core::chemical::ResidueType const & restype( residue_set->name_map( "ALA" ) );
        residue = core::conformation::ResidueFactory::create_residue( restype );

				//LSM taken from test structure
				VectorPair lsm_test(
						Vector(57.718,-16.524,29.748),
						Vector(0.16551,-0.93345,0.21081));

				LSMSearchPattern test_pattern(
						lsm_test,
						90,
						1,
						2,
						1,
						2);
				utility::vector1<RT> points_test = test_pattern.Searchpoints();

				foreach(RT transform, points_test)
				{
					core::pose::Pose testPose(targetPose);

					place_and_assert_transform(testPose, residue, transform);
				}
			}
  };
}
