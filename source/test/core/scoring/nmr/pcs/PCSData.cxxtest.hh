// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSData.cxxtest.hh
/// @brief   unit test for class PCSData that stores and handles all PCS data for all tagging sites and all lanthanides
/// @details Last modified: 07/19/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pcs/PCSData.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
#include <utility>

static basic::Tracer TR("core.scoring.nmr.pcs.PCSData.cxxtest");

class PCSDataTests : public CxxTest::TestSuite {

private:
	core::pose::Pose psR2_;
	core::pose::Pose cam_;

public:
	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(psR2_, "core/scoring/nmr/pcs/2ksy_chainA_renum.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		psR2_.clear();
		cam_.clear();
	}

	/// @brief Test the creation of PCSData
	void test_PCSData_instantiation() {
		using namespace core::scoring::nmr::pcs;

		PCSData pcs_data_psR2("core/scoring/nmr/pcs/pcs_data_2ksy_input_file.txt", psR2_);
		pcs_data_psR2.show(TR.Debug);
		TS_ASSERT_EQUALS(pcs_data_psR2.get_number_tags(), 5);

		PCSData pcs_data_cam("core/scoring/nmr/pcs/pcs_data_1cdl_mtsl_input_file.txt", cam_);
		pcs_data_cam.show(TR.Debug);
		TS_ASSERT_EQUALS(pcs_data_cam.get_number_tags(), 2);
	}

	/// @brief Test the creation of PCSData through the common interface of the NMRDataFactory.
	void test_PCSData_creation_by_NMRDataFactory() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		PCSDataOP pcs_data_psR2(utility::pointer::static_pointer_cast< PCSData >( NMRDataFactory::get_instance()->get_nmr_data("PCS", "core/scoring/nmr/pcs/pcs_data_2ksy_input_file.txt", psR2_) ) );
		pcs_data_psR2->show(TR.Debug);
		TS_ASSERT_EQUALS(pcs_data_psR2->get_number_tags(), 5);

		PCSDataOP pcs_data_cam(utility::pointer::static_pointer_cast< PCSData >( NMRDataFactory::get_instance()->get_nmr_data("PCS", "core/scoring/nmr/pcs/pcs_data_1cdl_mtsl_input_file.txt", cam_) ) );
		pcs_data_cam->show(TR.Debug);
		TS_ASSERT_EQUALS(pcs_data_cam->get_number_tags(), 2);
	}

	/// @brief Test PCS score calculation with grid search
	void test_PCSData_score_calculation_with_gridsearch() {
		using namespace core::scoring::nmr::pcs;

		TR << "Testing PCSData score calculation with grid search" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		PCSData pcs_data("core/scoring/nmr/pcs/pcs_data_2ksy_input_file.txt", psR2_);
		utility::vector1< utility::vector1< PCSTensorCOP > > tensors_all_tags_and_lanthanides(pcs_data.get_number_tags());
		utility::vector1< core::Real > scores_all_tags(pcs_data.get_number_tags());
		core::Real total_score = pcs_data.compute_score_all_tags(psR2_, scores_all_tags, tensors_all_tags_and_lanthanides);

		TR.Debug << "Total PCS score is " << total_score << std::endl;

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().z(),  -9.723, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().x(), 29.064, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().y(), 19.570, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().z(),  4.641, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().x(), 29.067, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().y(), 19.568, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().z(),  4.642, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().x(), 29.063, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().y(), 19.570, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().z(),  4.640, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().x(), 29.064, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().y(), 19.567, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().z(),  4.640, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().x(), -15.006, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().y(),   6.369, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().z(),  19.661, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().x(), -15.006, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().y(),   6.369, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().z(),  19.661, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().y(),   5.772, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().x(), -14.265, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().y(),   5.770, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().z(),  18.980, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().x(),  4.612, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().y(),  0.960, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().z(), 17.233, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().x(),  4.612, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().y(),  0.960, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().z(), 17.234, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().z(),  18.978, 1e-1);
	}

	/// @brief Test PCS score calculation with spinlabel
	void test_PCSData_score_calculation_with_spinlabel() {
		using namespace core::scoring::nmr::pcs;

		// Additional options, sets also fixed RG seed
		core_init_with_additional_options("-nmr:spinlabel:highres_conformer_filter_type DISTANCE");

		TR << "Testing PCSData score calculation with spinlabel" << std::endl;

		PCSData pcs_data("core/scoring/nmr/pcs/pcs_data_1cdl_mtsl_input_file.txt", cam_);
		utility::vector1< utility::vector1< PCSTensorCOP > > tensors_all_tags_and_lanthanides(pcs_data.get_number_tags());
		utility::vector1< core::Real > scores_all_tags(pcs_data.get_number_tags());
		core::Real total_score = pcs_data.compute_score_all_tags(cam_, scores_all_tags, tensors_all_tags_and_lanthanides);

		TR.Debug << "Total PCS score is " << total_score << std::endl;

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().x(), 85.9314, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().y(), 32.8311, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().z(),  9.3183, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().x(), 85.9314, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().y(), 32.8311, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().z(),  9.3183, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().x(), 85.9314, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().y(), 32.8311, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().z(),  9.3183, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().x(), 85.9314, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().y(), 32.8311, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().z(),  9.3183, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().x(), 78.3018, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().y(),  8.4364, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().z(), -4.1476, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().x(), 78.3018, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().y(),  8.4364, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().z(), -4.1476, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().x(), 78.3018, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().y(),  8.4364, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().z(), -4.1476, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().x(), 78.3018, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().y(),  8.4364, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().z(), -4.1476, 1e-1);

	}

	/// @brief Test PCS score calculation for scaled PCS data and grid search
	void test_scaled_PCSData_score_calculation() {
		using namespace core::scoring::nmr::pcs;

		// Additional options, sets also fixed RG seed
		core_init_with_additional_options("-nmr:pcs:normalize_data");

		TR << "Testing PCSData score calculation with scaled data and grid search" << std::endl;

		PCSData pcs_data("core/scoring/nmr/pcs/pcs_data_2ksy_input_file.txt", psR2_);
		utility::vector1< utility::vector1< PCSTensorCOP > > tensors_all_tags_and_lanthanides(pcs_data.get_number_tags());
		utility::vector1< core::Real > scores_all_tags(pcs_data.get_number_tags());
		core::Real total_score = pcs_data.compute_score_all_tags(psR2_, scores_all_tags, tensors_all_tags_and_lanthanides);

		TR.Debug << "Total PCS score is " << total_score << std::endl;

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().z(),  -9.723, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().x(),  22.496, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().y(), -19.107, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][4]->get_metal_center().z(),  -9.723, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().x(), 29.064, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().y(), 19.570, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][1]->get_metal_center().z(),  4.641, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().x(), 29.067, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().y(), 19.568, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][2]->get_metal_center().z(),  4.642, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().x(), 29.063, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().y(), 19.570, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][3]->get_metal_center().z(),  4.640, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().x(), 29.064, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().y(), 19.567, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[2][4]->get_metal_center().z(),  4.640, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().x(), -15.006, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().y(),   6.369, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][1]->get_metal_center().z(),  19.661, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().x(), -15.006, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().y(),   6.369, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][2]->get_metal_center().z(),  19.661, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().y(),   5.772, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][3]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().x(), -14.265, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().y(),   5.770, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[3][4]->get_metal_center().z(),  18.980, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().x(),  4.612, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().y(),  0.960, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][1]->get_metal_center().z(), 17.233, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().x(),  4.612, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().y(),  0.960, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[4][2]->get_metal_center().z(), 17.234, 1e-1);

		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][1]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][2]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][3]->get_metal_center().z(),  18.978, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().x(), -14.264, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().y(),   5.773, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[5][4]->get_metal_center().z(),  18.978, 1e-1);
	}
};
