// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   
/// @brief  test suite for protocols::toolbox::superimpose
/// @author Alex Ford (fordas@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package Headers
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <core/kinematics/Jump.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>
#include <protocols/toolbox/superimpose.hh>
#include <string>

static basic::Tracer TR("protocols.toolbox.superimpose.cxxtest");

class ToolboxSuperimposeTests : public CxxTest::TestSuite{


public: //typedef


	typedef core::conformation::Residue Residue;
	typedef core::pose::Pose Pose;

	typedef numeric::xyzVector<core::Real> Vector;
	typedef numeric::xyzMatrix<core::Real> Matrix;
	
public: //

	ToolboxSuperimposeTests() {};

  core::pose::PoseOP testpose_a;
  core::pose::PoseOP testpose_b;
  core::pose::PoseOP testpose_c;

  core::Real translation_magnitude;

	void setUp()
  {
		core_init_with_additional_options( "" );

    translation_magnitude = 10;

    // Generate test poses of 2 residues w/ N->N jumps
    // A translating 10A X
    // B translating 20A X
    // C translating 10A X, 10 A Y
    testpose_a = new core::pose::Pose();
    core::pose::make_pose_from_sequence(
      *testpose_a,
      "A",
      *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ));

    testpose_a->append_residue_by_jump( testpose_a->residue(1), 1, "N", "N", true);
    core::kinematics::Jump j_a = testpose_a->jump(1);
    j_a.translation_along_axis(
        testpose_a->conformation().upstream_jump_stub(1),
        Vector(1, 0, 0),
        translation_magnitude);
    testpose_a->set_jump_now(1, j_a);

    testpose_b = new core::pose::Pose();
    core::pose::make_pose_from_sequence(
      *testpose_b,
      "A",
      *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ));

    testpose_b->append_residue_by_jump( testpose_b->residue(1), 1, "N", "N", true);
    core::kinematics::Jump j_b = testpose_b->jump(1);
    j_b.translation_along_axis(
        testpose_b->conformation().upstream_jump_stub(1),
        Vector(1, 0, 0),
        2 * translation_magnitude);
    testpose_b->set_jump_now(1, j_b);

    testpose_c = new core::pose::Pose();
    core::pose::make_pose_from_sequence(
      *testpose_c,
      "A",
      *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ));

    testpose_c->append_residue_by_jump( testpose_c->residue(1), 1, "N", "N", true);
    core::kinematics::Jump j_c = testpose_c->jump(1);
    
    j_c.translation_along_axis(
        testpose_c->conformation().upstream_jump_stub(1),
        Vector(1, 0, 0),
        translation_magnitude);

    j_c.translation_along_axis(
        testpose_c->conformation().upstream_jump_stub(1),
        Vector(0, 1, 0),
        translation_magnitude);
    testpose_c->set_jump_now(1, j_c);
	}

  void test_superposition_by_jump()
  {
    using namespace protocols::toolbox;
    Vector ab_to_src, ab_to_fit;
    Matrix ab_rot;
		Vector result_transform;

		utility::vector1<Vector> testpose_a_all_coords = get_all_residue_bb_cords(*testpose_a);
		utility::vector1<Vector> testpose_b_all_coords = get_all_residue_bb_cords(*testpose_b);

    superposition_transform(
        testpose_a_all_coords, testpose_b_all_coords,
        ab_rot, ab_to_src, ab_to_fit);

		result_transform = ab_to_src - ab_to_fit;
    TS_ASSERT_DELTA(result_transform, Vector(translation_magnitude / 2, 0, 0), 1e-6);
    TS_ASSERT_DELTA(ab_rot, Matrix::identity(), 1e-6);

		core::pose::Pose testpose_a_cp1(*testpose_a);
		apply_superposition_transform_to_jump(testpose_a_cp1, 1, ab_rot, ab_to_src, ab_to_fit);

		TS_ASSERT_DELTA(
				testpose_a_cp1.residue(1).xyz("N"),
				testpose_a->residue(1).xyz("N"),
				1e-6);

		TS_ASSERT_DELTA(
				testpose_a_cp1.residue(2).xyz("N"),
				testpose_a->residue(2).xyz("N") + Vector(translation_magnitude / 2, 0, 0),
				1e-6);

		utility::vector1<Vector> testpose_a_2_coords = get_residue_2_bb_cords(*testpose_a);
		utility::vector1<Vector> testpose_b_2_coords = get_residue_2_bb_cords(*testpose_b);

    superposition_transform(
        testpose_a_2_coords,
        testpose_b_2_coords,
        ab_rot, ab_to_src, ab_to_fit);

		result_transform = ab_to_src - ab_to_fit;
    TS_ASSERT_DELTA(result_transform, Vector(translation_magnitude, 0, 0), 1e-6);
    TS_ASSERT_DELTA(ab_rot, Matrix::identity(), 1e-6);

		core::pose::Pose testpose_a_cp2(*testpose_a);
		apply_superposition_transform_to_jump(testpose_a_cp2, 1, ab_rot, ab_to_src, ab_to_fit);

		TS_ASSERT_DELTA(
				testpose_a_cp2.residue(1).xyz("N"),
				testpose_a->residue(1).xyz("N"),
				1e-6);

		TS_ASSERT_DELTA(
				testpose_a_cp2.residue(2).xyz("N"),
				testpose_a->residue(2).xyz("N") + Vector(translation_magnitude, 0, 0),
				1e-6);
  }

  void test_superposition()
  {
    using namespace protocols::toolbox;
    Vector ab_to_src, ab_to_fit;
    Matrix ab_rot;
		Vector result_transform;

		utility::vector1<Vector> testpose_a_all_coords = get_all_residue_bb_cords(*testpose_a);
		utility::vector1<Vector> testpose_b_all_coords = get_all_residue_bb_cords(*testpose_b);

    superposition_transform(
        testpose_a_all_coords, testpose_b_all_coords,
        ab_rot, ab_to_src, ab_to_fit);

		result_transform = ab_to_src - ab_to_fit;
    TS_ASSERT_DELTA(result_transform, Vector(translation_magnitude / 2, 0, 0), 1e-6);
    TS_ASSERT_DELTA(ab_rot, Matrix::identity(), 1e-6);

		core::pose::Pose testpose_a_cp1(*testpose_a);
		apply_superposition_transform(testpose_a_cp1, ab_rot, ab_to_src, ab_to_fit);

		TS_ASSERT_DELTA(
				testpose_a_cp1.residue(1).xyz("N"),
				testpose_a->residue(1).xyz("N") + Vector(translation_magnitude / 2, 0, 0),
				1e-6);

		TS_ASSERT_DELTA(
				testpose_a_cp1.residue(2).xyz("N"),
				testpose_a->residue(2).xyz("N") + Vector(translation_magnitude / 2, 0, 0),
				1e-6);

		utility::vector1<Vector> testpose_a_2_coords = get_residue_2_bb_cords(*testpose_a);
		utility::vector1<Vector> testpose_b_2_coords = get_residue_2_bb_cords(*testpose_b);

    superposition_transform(
        testpose_a_2_coords,
        testpose_b_2_coords,
        ab_rot, ab_to_src, ab_to_fit);

		result_transform = ab_to_src - ab_to_fit;
    TS_ASSERT_DELTA(result_transform, Vector(translation_magnitude, 0, 0), 1e-6);
    TS_ASSERT_DELTA(ab_rot, Matrix::identity(), 1e-6);

		core::pose::Pose testpose_a_cp2(*testpose_a);
		apply_superposition_transform(testpose_a_cp2, ab_rot, ab_to_src, ab_to_fit);

		TS_ASSERT_DELTA(
				testpose_a_cp2.residue(1).xyz("N"),
				testpose_a->residue(1).xyz("N") + Vector(translation_magnitude, 0, 0),
				1e-6);

		TS_ASSERT_DELTA(
				testpose_a_cp2.residue(2).xyz("N"),
				testpose_a->residue(2).xyz("N") + Vector(translation_magnitude, 0, 0),
				1e-6);
  }

  utility::vector1<Vector> get_all_residue_bb_cords(core::pose::Pose & pose)
  {
    utility::vector1<Vector> result;
    result.push_back(pose.residue(1).xyz("N"));
    result.push_back(pose.residue(1).xyz("C"));
    result.push_back(pose.residue(1).xyz("CA"));
    result.push_back(pose.residue(2).xyz("N"));
    result.push_back(pose.residue(2).xyz("C"));
    result.push_back(pose.residue(2).xyz("CA"));
    return result;
  }

  utility::vector1<Vector> get_residue_2_bb_cords(core::pose::Pose & pose)
  {
    utility::vector1<Vector> result;
    result.push_back(pose.residue(2).xyz("N"));
    result.push_back(pose.residue(2).xyz("C"));
    result.push_back(pose.residue(2).xyz("CA"));
    return result;
  }

  core::pose::PoseOP test_pose_a;
  core::pose::PoseOP test_pose_b;
};

