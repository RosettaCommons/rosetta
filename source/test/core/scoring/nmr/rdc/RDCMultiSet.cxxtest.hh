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

static basic::Tracer TR("core.scoring.nmr.rdc.RDCMultiSet.cxxtest");

class RDCMultiSetTests : public CxxTest::TestSuite {

private: // Data

	core::pose::Pose g_;
	utility::vector1< std::string > bicelle_exp_input_;
	utility::vector1< std::string > tmv107_exp_input_;

public: // Test functions

	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		// RDC data are not pre-scaled yet and we have to correct the sign of the 15N gyromagnetic ratio
		core_init_with_additional_options("-nmr:rdc:correct_sign true -nmr:rdc:normalization_type NH");

		// Load pose from pdb
		core::import_pose::pose_from_file(g_, "core/scoring/nmr/rdc/gb1.pdb", core::import_pose::PDB_file);

		bicelle_exp_input_.push_back("core/scoring/nmr/rdc/bicelles_nh.txt");
		bicelle_exp_input_.push_back("core/scoring/nmr/rdc/bicelles_nc.txt");
		bicelle_exp_input_.push_back("core/scoring/nmr/rdc/bicelles_hnc.txt");

		tmv107_exp_input_.push_back("core/scoring/nmr/rdc/tmv107_nh.txt");
		tmv107_exp_input_.push_back("core/scoring/nmr/rdc/tmv107_nc.txt");
		tmv107_exp_input_.push_back("core/scoring/nmr/rdc/tmv107_hnc.txt");
	}

	void tearDown() {
		g_.clear();
		bicelle_exp_input_.clear();
		tmv107_exp_input_.clear();
	}

	/// @brief test RDCMultiSet creation and getters
	void test_RDCMultiSet_instantiation() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCMultiSet instantiation" << std::endl;

		RDCMultiSet multiset_1(bicelle_exp_input_, "Bicelles", g_);
		RDCMultiSet multiset_2(tmv107_exp_input_, "TMV107", g_, 1.0, "NLSDA");

		// Get alignment label
		TS_ASSERT_EQUALS(multiset_1.get_alignment_medium_label(), "Bicelles");
		TS_ASSERT_EQUALS(multiset_2.get_alignment_medium_label(), "TMV107");

		// Get number of experiments
		TS_ASSERT_EQUALS(multiset_1.get_number_experiments(), 3);
		TS_ASSERT_EQUALS(multiset_2.get_number_experiments(), 3);

		// Get total number of RDCs
		TS_ASSERT_EQUALS(multiset_1.get_total_number_rdc(), 147);
		TS_ASSERT_EQUALS(multiset_2.get_total_number_rdc(), 152);

		// Get number of RDCs for each experiment
		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[1]->get_number_rdc(), 46);
		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[2]->get_number_rdc(), 50);
		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[3]->get_number_rdc(), 51);
		TS_ASSERT_EQUALS(multiset_2.get_rdc_singleset_vec()[1]->get_number_rdc(), 48);
		TS_ASSERT_EQUALS(multiset_2.get_rdc_singleset_vec()[2]->get_number_rdc(), 51);
		TS_ASSERT_EQUALS(multiset_2.get_rdc_singleset_vec()[3]->get_number_rdc(), 53);

		// Test accessing the RDCSingleSets
		multiset_1.get_rdc_singleset_vec()[1]->set_single_rdc_weighting_scheme("CONST");
		multiset_1.get_rdc_singleset_vec()[2]->set_single_rdc_weighting_scheme("SIGMA");
		multiset_1.get_rdc_singleset_vec()[3]->set_single_rdc_weighting_scheme("OBSIG");
		multiset_2.get_rdc_singleset_vec()[1]->set_weight(1.5);
		multiset_2.get_rdc_singleset_vec()[2]->set_weight(0.5);
		multiset_2.get_rdc_singleset_vec()[3]->set_weight(2.5);

		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[1]->get_single_rdc_weighting_scheme(), CONST);
		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[2]->get_single_rdc_weighting_scheme(), SIGMA);
		TS_ASSERT_EQUALS(multiset_1.get_rdc_singleset_vec()[3]->get_single_rdc_weighting_scheme(), OBSIG);
		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[1]->get_weight(), 1.5, 1.0e-6);
		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[2]->get_weight(), 0.5, 1.0e-6);
		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[3]->get_weight(), 2.5, 1.0e-6);

		// Test accessing single RDC values
		// These should be normalized now
		core::Size offset1(multiset_1.get_rdc_singleset_vec()[1]->get_number_rdc());
		core::Size offset2 = offset1 + multiset_1.get_rdc_singleset_vec()[2]->get_number_rdc();
		TS_ASSERT_DELTA(multiset_1.get_rdc_values()(1), 2.269, 1e-3);
		TS_ASSERT_DELTA(multiset_1.get_rdc_values()(offset1+2), -6.4217,1e-3);
		TS_ASSERT_DELTA(multiset_1.get_rdc_values()(offset2+3), -3.8474,1e-3);

		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[1]->get_single_rdc_vec()[1].get_rdc_exp(), -3.733, 1e-6);
		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[2]->get_single_rdc_vec()[2].get_rdc_exp(),  0.19033432, 1e-6);
		TS_ASSERT_DELTA(multiset_2.get_rdc_singleset_vec()[3]->get_single_rdc_vec()[3].get_rdc_exp(),  7.47760395, 1e-6);

	}

	/// @brief test RDC tensor and score calculation with SVD
	void test_score_calculation_by_svd() {
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCMultiSet score calculation with SVD" << std::endl;

		// Create RDCMultiSets, set singleset weights and weighting scheme
		// Update single rdc weights, spin coordinates and matrix A before tensor and score calculation
		RDCMultiSet multiset_1(bicelle_exp_input_, "Bicelles", g_, 1.0, "SVD");
		for ( core::Size i = 1; i <= multiset_1.get_number_experiments(); ++i ) {
			multiset_1.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_1.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_1.update_single_rdc_weighting();
		RDCMultiSet multiset_2(tmv107_exp_input_, "TMV107", g_, 1.0, "SVD");
		for ( core::Size i = 1; i <= multiset_2.get_number_experiments(); ++i ) {
			multiset_2.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_2.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_2.update_single_rdc_weighting();

		multiset_1.update_spin_coordinates(g_);
		multiset_1.update_matrix_A();
		multiset_2.update_spin_coordinates(g_);
		multiset_2.update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score1 = multiset_1.solve_tensor_and_compute_score_by_svd();
		core::Real score2 = multiset_2.solve_tensor_and_compute_score_by_svd();

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_1.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score1 << std::endl;
		RDCTensorOP tensor1(multiset_1.get_tensor());
		tensor1->show_tensor_stats(TR.Debug);
		tensor1->set_Dmax(21584.630);
		tensor1->diagonalize_tensor();
		tensor1->reorder_tensor();
		TS_ASSERT_DELTA(tensor1->get_T_xx(),  -5.128E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_xy(),  -1.008E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_xz(),   1.387E-05, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_yy(),   1.448E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_yz(),   5.319E-05, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_ax(),    -7.927E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_rh(),    -2.318E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_Da(),    -8.555, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_R(),      0.292, 1e-1);
		tensor1->show_tensor_stats(TR.Debug);

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_2.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score2 << std::endl;
		RDCTensorOP tensor2(multiset_2.get_tensor());
		tensor2->show_tensor_stats(TR.Debug);
		tensor2->set_Dmax(21584.630);
		tensor2->diagonalize_tensor();
		tensor2->reorder_tensor();
		TS_ASSERT_DELTA(tensor2->get_T_xx(),   3.166E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_T_xy(),  -4.489E-05, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_T_xz(),   7.691E-05, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_T_yy(),  -2.289E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_T_yz(),  -1.683E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_ax(),     5.123E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_rh(),     3.401E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_Da(),     5.528, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_R(),      0.664, 1e-1);
		tensor2->show_tensor_stats(TR.Debug);
	}

	/// @brief test RDC tensor and score calculation with NLS
	void test_score_calculation_by_nls() {
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCMultiSet score calculation with NLS fitting" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// Create RDCMultiSets, set singleset weights and weighting scheme
		// Update single rdc weights, spin coordinates and matrix A before tensor and score calculation
		RDCMultiSet multiset_1(bicelle_exp_input_, "Bicelles", g_, 1.0, "NLS");
		for ( core::Size i = 1; i <= multiset_1.get_number_experiments(); ++i ) {
			multiset_1.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_1.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_1.update_single_rdc_weighting();
		RDCMultiSet multiset_2(tmv107_exp_input_, "TMV107", g_, 1.0, "NLS");
		for ( core::Size i = 1; i <= multiset_2.get_number_experiments(); ++i ) {
			multiset_2.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_2.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_2.update_single_rdc_weighting();

		multiset_1.update_spin_coordinates(g_);
		multiset_1.update_matrix_A();
		multiset_2.update_spin_coordinates(g_);
		multiset_2.update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score1 = multiset_1.solve_tensor_and_compute_score_by_nls();
		core::Real score2 = multiset_2.solve_tensor_and_compute_score_by_nls();

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_1.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score1 << std::endl;
		RDCTensorOP tensor1(multiset_1.get_tensor());
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		//Unconstrained nls fitting of this particular RDC dataset
		//doesn't lead to the best possible solution as opposed to
		//using a fixed value for Da or R (see test functions below).
		//This is why I decided not to test it here
		//TS_ASSERT_DELTA(tensor1->get_ax(),    -7.927E-04, 1e-5);
		//TS_ASSERT_DELTA(tensor1->get_rh(),    -2.318E-04, 1e-5);
		//TS_ASSERT_DELTA(tensor1->get_Da(),    -8.555, 1e-1);
		//TS_ASSERT_DELTA(tensor1->get_R(),      0.292, 1e-1);

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_2.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score2 << std::endl;
		RDCTensorOP tensor2(multiset_2.get_tensor());
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		//Unconstrained nls fitting of this particular RDC dataset
		//doesn't lead to the best possible solution as opposed to
		//using a fixed value for Da or R (see test functions below).
		//This is why I decided not to test it here
		//TS_ASSERT_DELTA(tensor2->get_ax(),     5.123E-04, 1e-5);
		//TS_ASSERT_DELTA(tensor2->get_rh(),     3.401E-04, 1e-5);
		//TS_ASSERT_DELTA(tensor2->get_Da(),     5.528, 1e-1);
		//TS_ASSERT_DELTA(tensor2->get_R(),      0.664, 1e-1);
	}

	/// @brief test RDC tensor and score calculation with NLS and predefined value for Da
	void test_score_calculation_by_nls_da() {
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCMultiSet score calculation with NLS fitting and fixed Da" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// Create RDCMultiSets, set singleset weights and weighting scheme
		// Update single rdc weights, spin coordinates and matrix A before tensor and score calculation
		RDCMultiSet multiset_1(bicelle_exp_input_, "Bicelles", g_, 1.0, "NLSDA");
		for ( core::Size i = 1; i <= multiset_1.get_number_experiments(); ++i ) {
			multiset_1.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_1.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_1.update_single_rdc_weighting();
		utility::vector1<core::Real> fixed_tensor_values(6);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360.0;     // alpha, beta, gamma
		fixed_tensor_values[4] = 99.9; fixed_tensor_values[5] = 0.6; fixed_tensor_values[6] = 21584.63; // Da, R, Dmax

		fixed_tensor_values[4] = -8.555;
		RDCTensorOP tensor1( new RDCTensor );
		tensor1->set_tensor_in_pas(fixed_tensor_values);
		multiset_1.set_tensor(tensor1);

		RDCMultiSet multiset_2(tmv107_exp_input_, "TMV107", g_, 1.0, "NLSDA");
		for ( core::Size i = 1; i <= multiset_2.get_number_experiments(); ++i ) {
			multiset_2.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_2.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_2.update_single_rdc_weighting();
		fixed_tensor_values[4] = 5.528;
		RDCTensorOP tensor2( new RDCTensor );
		tensor2->set_tensor_in_pas(fixed_tensor_values);
		multiset_2.set_tensor(tensor2);

		multiset_1.update_spin_coordinates(g_);
		multiset_1.update_matrix_A();
		multiset_2.update_spin_coordinates(g_);
		multiset_2.update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score1 = multiset_1.solve_tensor_and_compute_score_by_nls();
		core::Real score2 = multiset_2.solve_tensor_and_compute_score_by_nls();

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_1.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score1 << std::endl;
		tensor1 = multiset_1.get_tensor();
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor1->get_ax(),     -7.927E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_rh(),     -2.318E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_Da(),     -8.555, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_R(),       0.292, 1e-1);

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_2.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score2 << std::endl;
		tensor2 = multiset_2.get_tensor();
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor2->get_ax(),    5.122E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_rh(),    3.401E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_Da(),    5.528, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_R(),     0.664, 1e-1);
	}

	/// @brief test RDC tensor and score calculation with NLS and predefined value for R
	void test_score_calculation_by_nls_r() {
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCMultiSet score calculation with NLS fitting and fixed R" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// Create RDCMultiSets, set singleset weights and weighting scheme
		// Update single rdc weights, spin coordinates and matrix A before tensor and score calculation
		RDCMultiSet multiset_1(bicelle_exp_input_, "Bicelles", g_, 1.0, "NLSR");
		for ( core::Size i = 1; i <= multiset_1.get_number_experiments(); ++i ) {
			multiset_1.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_1.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_1.update_single_rdc_weighting();
		utility::vector1<core::Real> fixed_tensor_values(6);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360.0;     // alpha, beta, gamma
		fixed_tensor_values[4] = 99.9; fixed_tensor_values[5] = 0.6; fixed_tensor_values[6] = 21584.63; // Da, R, Dmax

		fixed_tensor_values[5] = 0.292;
		RDCTensorOP tensor1( new RDCTensor );
		tensor1->set_tensor_in_pas(fixed_tensor_values);
		multiset_1.set_tensor(tensor1);

		RDCMultiSet multiset_2(tmv107_exp_input_, "TMV107", g_, 1.0, "NLSR");
		for ( core::Size i = 1; i <= multiset_2.get_number_experiments(); ++i ) {
			multiset_2.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_2.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_2.update_single_rdc_weighting();
		fixed_tensor_values[5] = 0.664;
		RDCTensorOP tensor2( new RDCTensor );
		tensor2->set_tensor_in_pas(fixed_tensor_values);
		multiset_2.set_tensor(tensor2);

		multiset_1.update_spin_coordinates(g_);
		multiset_1.update_matrix_A();
		multiset_2.update_spin_coordinates(g_);
		multiset_2.update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score1 = multiset_1.solve_tensor_and_compute_score_by_nls();
		core::Real score2 = multiset_2.solve_tensor_and_compute_score_by_nls();

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_1.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score1 << std::endl;
		tensor1 = multiset_1.get_tensor();
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor1->get_ax(),      -7.928E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_rh(),      -2.315E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_Da(),      -8.556, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_R(),        0.292, 1e-1);

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_2.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score2 << std::endl;
		tensor2 = multiset_2.get_tensor();
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor2->get_ax(),    -5.120E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_rh(),    -3.400E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_Da(),    -5.526, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_R(),      0.664, 1e-1);
	}

	/// @brief test RDC tensor and score calculation with pre-normalized RDC data.
	///        thus, we turn the normalization off. Also, we don't have to correct
	///        the sign of the 15N gyromagnetic ratio because the sign of the
	///        NCO and COHN couplings is already flipped.
	void test_score_calculation_with_normalized_data() {
		using namespace core::scoring::nmr::rdc;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:normalization_type none -nmr:rdc:correct_sign false");

		TR << "Testing RDCMultiSet score calculation for normalized data" << std::endl;

		utility::vector1< std::string > bicelle_normalized_data;
		bicelle_normalized_data.push_back("core/scoring/nmr/rdc/bicelles_nh_norm_toNH.txt");
		bicelle_normalized_data.push_back("core/scoring/nmr/rdc/bicelles_nc_norm_toNH.txt");
		bicelle_normalized_data.push_back("core/scoring/nmr/rdc/bicelles_hnc_norm_toNH.txt");

		utility::vector1< std::string > tmv107_normalized_data;
		tmv107_normalized_data.push_back("core/scoring/nmr/rdc/tmv107_nh_norm_toNH.txt");
		tmv107_normalized_data.push_back("core/scoring/nmr/rdc/tmv107_nc_norm_toNH.txt");
		tmv107_normalized_data.push_back("core/scoring/nmr/rdc/tmv107_hnc_norm_toNH.txt");

		// Create RDCMultiSets, set singleset weights and weighting scheme
		// Update single rdc weights, spin coordinates and matrix A before tensor and score calculation
		RDCMultiSet multiset_1(bicelle_normalized_data, "Bicelles", g_, 1.0, "SVD");
		for ( core::Size i = 1; i <= multiset_1.get_number_experiments(); ++i ) {
			multiset_1.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_1.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_1.update_single_rdc_weighting();
		RDCMultiSet multiset_2(tmv107_normalized_data, "TMV107", g_, 1.0, "NLSDA");
		for ( core::Size i = 1; i <= multiset_2.get_number_experiments(); ++i ) {
			multiset_2.get_rdc_singleset_vec()[i]->set_single_rdc_weighting_scheme("CONST");
			multiset_2.get_rdc_singleset_vec()[i]->set_weight(1.0);
		}
		multiset_2.update_single_rdc_weighting();
		utility::vector1<core::Real> fixed_tensor_values(6);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360.0;      // alpha, beta, gamma
		fixed_tensor_values[4] = 5.528; fixed_tensor_values[5] = 0.6; fixed_tensor_values[6] = 21584.63; // Da, R, Dmax
		RDCTensorOP tensor2( new RDCTensor );
		tensor2->set_tensor_in_pas(fixed_tensor_values);
		multiset_2.set_tensor(tensor2);

		multiset_1.update_spin_coordinates(g_);
		multiset_1.update_matrix_A();
		multiset_2.update_spin_coordinates(g_);
		multiset_2.update_matrix_A();

		// Calculate and show RDC scores and tensors
		core::Real score1 = multiset_1.solve_tensor_and_compute_score_by_svd();
		core::Real score2 = multiset_2.solve_tensor_and_compute_score_by_nls();

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_1.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score1 << std::endl;
		RDCTensorOP tensor1(multiset_1.get_tensor());
		tensor1->show_tensor_stats(TR.Debug);
		tensor1->set_Dmax(21584.630);
		tensor1->diagonalize_tensor();
		tensor1->reorder_tensor();
		tensor1->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor1->get_T_xx(),  -5.128E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_xy(),  -1.009E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_xz(),   1.387E-05, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_yy(),   1.448E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_T_yz(),   5.319E-05, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_ax(),    -7.927E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_rh(),    -2.318E-04, 1e-5);
		TS_ASSERT_DELTA(tensor1->get_Da(),    -8.555, 1e-1);
		TS_ASSERT_DELTA(tensor1->get_R(),      0.292, 1e-1);

		TR.Debug << "Calculated score and alignment tensor for medium: " << multiset_2.get_alignment_medium_label() << std::endl;
		TR.Debug << "RDC score = " << score2 << std::endl;
		tensor2 = multiset_2.get_tensor();
		tensor2->reorder_tensor();
		tensor2->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor2->get_ax(),     5.122E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_rh(),     3.401E-04, 1e-5);
		TS_ASSERT_DELTA(tensor2->get_Da(),     5.528, 1e-1);
		TS_ASSERT_DELTA(tensor2->get_R(),      0.664, 1e-1);

	}
};
