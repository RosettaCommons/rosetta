// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCMultiSet_Symmetry.cxxtest.hh
/// @brief   unit test for class RDCMultiSet that stores and handles data of one RDC alignment medium
///          that can include experiments collected for different spin types
/// @details Last modified: 07/09/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/scoring/nmr/rdc/RDCTensor.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

static basic::Tracer TR("core.scoring.nmr.rdc.RDCMultiSet_Symmetry.cxxtest");

class RDCMultiSetSymmetryTests : public CxxTest::TestSuite {

private: // Data

	core::pose::Pose tolR_;
	utility::vector1< std::string > rdc_input_asu_;
	utility::vector1< std::string > rdc_input_all_su_;

public: // Test functions

	/// @brief Setup Test
	void setUp() {
		using namespace core::conformation::symmetry;
		using namespace core::pose;

		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(tolR_, "core/scoring/nmr/rdc/2jwk_A_model1_renamed_INPUT.pdb", core::import_pose::PDB_file);

		// Make symmetric pose
		SymmDataOP symmdata = SymmDataOP( new SymmData );
		symmdata->read_symmetry_data_from_file( "core/scoring/nmr/rdc/2jwk_AB_model1.symm" );
		core::pose::symmetry::make_symmetric_pose( tolR_, *symmdata );

		rdc_input_asu_.push_back("core/scoring/nmr/rdc/2jwk_caco_asu.txt");
		rdc_input_asu_.push_back("core/scoring/nmr/rdc/2jwk_caha_asu.txt");
		rdc_input_asu_.push_back("core/scoring/nmr/rdc/2jwk_nco_asu.txt");
		rdc_input_asu_.push_back("core/scoring/nmr/rdc/2jwk_nh_asu.txt");

		rdc_input_all_su_.push_back("core/scoring/nmr/rdc/2jwk_caco_all_su.txt");
		rdc_input_all_su_.push_back("core/scoring/nmr/rdc/2jwk_caha_all_su.txt");
		rdc_input_all_su_.push_back("core/scoring/nmr/rdc/2jwk_nco_all_su.txt");
		rdc_input_all_su_.push_back("core/scoring/nmr/rdc/2jwk_nh_all_su.txt");
	}

	void tearDown() {
		tolR_.clear();
		rdc_input_asu_.clear();
		rdc_input_all_su_.clear();
	}

	/// @brief Test RDC tensor and score calculation with SVD on a symmetric protein.
	///        In addition, residues in the pdb and rdc data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        The symmetric spins get deduced automatically from the symmetry info object
	///        and only the spins of the ASU are provided in the input file.
	void test_score_calculation_by_svd_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr::rdc;

		// Important: use this option to turn automatic symmetry handling on
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH -nmr:rdc:use_symmetry_calc");

		TR << "Testing RDCMultiSet score calculation for symmetric pose with SVD" << std::endl;

		// Setup rdc data
		RDCMultiSetOP rdc_data = RDCMultiSetOP( new RDCMultiSet(rdc_input_asu_, "phage_pf1", tolR_, 1.0, "SVD") );
		for ( core::Size i = 1; i <= rdc_data->get_number_experiments(); ++i ) {
			rdc_data->get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			rdc_data->get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		rdc_data->update_single_rdc_weighting();
		rdc_data->update_spin_coordinates(tolR_);
		rdc_data->update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score = rdc_data->solve_tensor_and_compute_score_by_svd();
		TR.Debug << "Calculated score and alignment tensor for medium: " << rdc_data->get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score << std::endl;
		RDCTensorOP tensor(rdc_data->get_tensor());
		tensor->show_tensor_stats(TR.Debug);
		tensor->set_Dmax(21584.630);
		tensor->diagonalize_tensor();
		tensor->reorder_tensor();
		tensor->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor->get_T_xx(),  -2.208E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_xy(),   1.922E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_xz(),   3.523E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_yy(),   2.179E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_yz(),   1.422E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_ax(),    -7.180E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_rh(),    -6.250E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_Da(),    -7.749, 1e-1);
		TS_ASSERT_DELTA(tensor->get_R(),      0.087, 1e-1);
	}

	/// @brief Test RDC tensor and score calculation with SVD on a symmetric protein.
	///        In addition, residues in the pdb and rdc data files do not start with 1,
	///        to test that conversion from pdb to pose numbering works properly.
	///        Here, the symmetric spins are explicitly declared in the input file
	///        and an identical set of rdcs is assigned the symmetric subunit.
	void test_score_calculation_by_svd_manual_symmetry_handling() {
		using namespace core::scoring::nmr::rdc;

		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH");

		TR << "Testing RDCMultiSet score calculation for symmetric pose with SVD" << std::endl;

		// Setup rdc data
		RDCMultiSetOP rdc_data = RDCMultiSetOP( new RDCMultiSet(rdc_input_all_su_, "phage_pf1", tolR_, 1.0, "SVD") );
		for ( core::Size i = 1; i <= rdc_data->get_number_experiments(); ++i ) {
			rdc_data->get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			rdc_data->get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		rdc_data->update_single_rdc_weighting();
		rdc_data->update_spin_coordinates(tolR_);
		rdc_data->update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score = rdc_data->solve_tensor_and_compute_score_by_svd();
		TR.Debug << "Calculated score and alignment tensor for medium: " << rdc_data->get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score << std::endl;
		RDCTensorOP tensor(rdc_data->get_tensor());
		tensor->show_tensor_stats(TR.Debug);
		tensor->set_Dmax(21584.630);
		tensor->diagonalize_tensor();
		tensor->reorder_tensor();
		tensor->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor->get_T_xx(),  -2.208E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_xy(),   1.922E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_xz(),   3.523E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_yy(),   2.179E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_T_yz(),   1.422E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_ax(),    -7.180E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_rh(),    -6.250E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_Da(),    -7.749, 1e-1);
		TS_ASSERT_DELTA(tensor->get_R(),      0.087, 1e-1);
	}

	/// @brief Same test as above but instead of SVD we use NLS for tensor
	///        and score determination.
	void test_score_calculation_by_nls_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr::rdc;

		// Important: use this option to turn automatic symmetry handling on
		// Sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH -nmr:rdc:use_symmetry_calc");

		TR << "Testing RDCMultiSet score calculation for symmetric pose with NLSDAR" << std::endl;

		// Setup rdc data
		RDCMultiSetOP rdc_data = RDCMultiSetOP( new RDCMultiSet(rdc_input_asu_, "phage_pf1", tolR_, 1.0, "NLSDAR") );
		for ( core::Size i = 1; i <= rdc_data->get_number_experiments(); ++i ) {
			rdc_data->get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("OBSIG");
			rdc_data->get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		rdc_data->update_single_rdc_weighting();
		utility::vector1<core::Real> fixed_tensor_values(6);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360.0;
		fixed_tensor_values[4] = -7.760; fixed_tensor_values[5] = 0.1; fixed_tensor_values[6] = 21584.63;
		RDCTensorOP tensor( new RDCTensor );
		tensor->set_tensor_in_pas(fixed_tensor_values);
		rdc_data->set_tensor(tensor);
		rdc_data->update_spin_coordinates(tolR_);
		rdc_data->update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score = rdc_data->solve_tensor_and_compute_score_by_nls();
		TR.Debug << "Calculated score and alignment tensor for medium: " << rdc_data->get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score << std::endl;
		tensor=rdc_data->get_tensor();
		tensor->reorder_tensor();
		tensor->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor->get_ax(),     -7.190E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_rh(),     -7.213E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_Da(),     -7.760, 1e-1);
		TS_ASSERT_DELTA(tensor->get_R(),       0.100, 1e-1);
	}

	/// @brief Same test as above but instead of SVD we use NLS for tensor
	///        and score determination.
	void test_score_calculation_by_nls_manual_symmetry_handling() {
		using namespace core::scoring::nmr::rdc;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH");

		TR << "Testing RDCMultiSet score calculation for symmetric pose with NLSDAR" << std::endl;

		// Setup rdc data
		RDCMultiSetOP rdc_data = RDCMultiSetOP( new RDCMultiSet(rdc_input_all_su_, "phage_pf1", tolR_, 1.0, "NLS") );
		for ( core::Size i = 1; i <= rdc_data->get_number_experiments(); ++i ) {
			rdc_data->get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("OBSIG");
			rdc_data->get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		rdc_data->update_single_rdc_weighting();
		utility::vector1<core::Real> fixed_tensor_values(6);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360.0;
		fixed_tensor_values[4] = -7.760; fixed_tensor_values[5] = 0.1; fixed_tensor_values[6] = 21584.63;
		RDCTensorOP tensor( new RDCTensor );
		tensor->set_tensor_in_pas(fixed_tensor_values);
		rdc_data->set_tensor(tensor);
		rdc_data->update_spin_coordinates(tolR_);
		rdc_data->update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score = rdc_data->solve_tensor_and_compute_score_by_nls();
		TR.Debug << "Calculated score and alignment tensor for medium: " << rdc_data->get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score << std::endl;
		tensor=rdc_data->get_tensor();
		tensor->reorder_tensor();
		tensor->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor->get_ax(),     -7.190E-04, 1e-5);
		TS_ASSERT_DELTA(tensor->get_rh(),     -7.213E-05, 1e-5);
		TS_ASSERT_DELTA(tensor->get_Da(),     -7.760, 1e-1);
		TS_ASSERT_DELTA(tensor->get_R(),       0.100, 1e-1);
	}
};
