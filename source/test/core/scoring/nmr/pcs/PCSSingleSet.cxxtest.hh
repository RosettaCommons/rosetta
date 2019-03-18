// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingleSet.cxxtest.hh
/// @brief   unit test for class PCSSingleSet that stores and handles data of one single PCS dataset (i.e. of one lanthanide ion)
/// @details Last modified: 07/09/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

static basic::Tracer TR("core.scoring.nmr.pcs.PCSSingleSet.cxxtest");

class PCSSingleSetTests : public CxxTest::TestSuite {

private:
	core::pose::Pose sh2_;

public:

	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(sh2_, "core/scoring/nmr/pcs/1x0n_model1_renum.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		sh2_.clear();
	}

	/// @brief test construction of PCSSingleSet object from PCS data file
	void test_PCSSingleSet_instantiation_and_getters() {
		using namespace core::scoring::nmr::pcs;

		// ctor with default values for single pcs weighting scheme and computation type
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_);

		// ctor with filename, weight, single pcs weighting scheme and computation type
		PCSSingleSet single_dataset_dy("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, 1.0, "CONST", "SVD");
		PCSSingleSet single_dataset_tm("core/scoring/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_, 1.0, "SIGMA", "NLS");
		PCSSingleSet single_dataset_er("core/scoring/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_, 1.0, "OBSIG", "NLSAXRH");

		// get number of pcs
		TS_ASSERT_EQUALS(single_dataset_tb.get_number_pcs(), 84);
		TS_ASSERT_EQUALS(single_dataset_dy.get_number_pcs(), 89);
		TS_ASSERT_EQUALS(single_dataset_tm.get_number_pcs(), 86);
		TS_ASSERT_EQUALS(single_dataset_er.get_number_pcs(), 94);

		// access pcs data
		TS_ASSERT_DELTA(single_dataset_tb.get_pcs_values()(1), -1.659, 1e-6);
		TS_ASSERT_DELTA(single_dataset_dy.get_pcs_values()(1),  0.001, 1e-6);
		TS_ASSERT_DELTA(single_dataset_tm.get_pcs_values()(1),  1.248, 1e-6);
		TS_ASSERT_DELTA(single_dataset_er.get_pcs_values()(1),  0.233, 1e-6);

		// get metal ion label
		TS_ASSERT_EQUALS(single_dataset_tb.get_metal_ion_label(), "Tb");
		TS_ASSERT_EQUALS(single_dataset_dy.get_metal_ion_label(), "Dy");
		TS_ASSERT_EQUALS(single_dataset_tm.get_metal_ion_label(), "Tm");
		TS_ASSERT_EQUALS(single_dataset_er.get_metal_ion_label(), "Er");

		// get computation type
		TS_ASSERT_EQUALS(single_dataset_tb.get_computation_type(), PCSSingleSet::SVD);
		TS_ASSERT_EQUALS(single_dataset_dy.get_computation_type(), PCSSingleSet::SVD);
		TS_ASSERT_EQUALS(single_dataset_tm.get_computation_type(), PCSSingleSet::NLS);
		TS_ASSERT_EQUALS(single_dataset_er.get_computation_type(), PCSSingleSet::NLSAXRH);

		// get pcs single weights
		TS_ASSERT_DELTA(single_dataset_tb.get_pcs_single_weights()(1),   1.0, 1e-6);
		TS_ASSERT_DELTA(single_dataset_dy.get_pcs_single_weights()(1),   1.0, 1e-6);
		TS_ASSERT_DELTA(single_dataset_tm.get_pcs_single_weights()(1), 100.0, 1e-6);
		TS_ASSERT_DELTA(single_dataset_er.get_pcs_single_weights()(1),  12.732240437, 1e-6);
	}

	/// @brief test construction of PCSSingleSet object from PCS data file and scaling of input pcs data
	void test_PCSSingleSet_instantiation_and_scaling() {
		using namespace core::scoring::nmr::pcs;
		core_init_with_additional_options("-nmr:pcs:normalize_data");

		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, 1.0, "CONST", "SVD");
		PCSSingleSet single_dataset_dy("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, 1.0, "SIGMA", "NLS");

		// normalized data
		TS_ASSERT(single_dataset_tb.normalized_data());
		TS_ASSERT(single_dataset_dy.normalized_data());

		// access scaled pcs values and weights
		TS_ASSERT_DELTA(single_dataset_tb.get_pcs_values()(1),         -3.3664957, 1e-6);
		TS_ASSERT_DELTA(single_dataset_tb.get_pcs_single_weights()(1),  1.0000000, 1e-6);
		TS_ASSERT_DELTA(single_dataset_dy.get_pcs_values()(1),          0.0020981, 1e-6);
		TS_ASSERT_DELTA(single_dataset_dy.get_pcs_single_weights()(1), 22.7157012, 1e-6);
	}

	/// @brief test PCS tensor calculation with grid search and SVD
	void test_calc_tensor_and_score_svd() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;

		TR << "Testing PCS tensor and score calculation by SVD" << std::endl;

		// create pcs dataset
		core::Real weight(1.0);
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, weight, "CONST", "SVD");

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearch grid_searcher(grid_atom1, grid_atom2, sh2_, 15, 1, 0, 25);

		// set grid search center coordinates
		grid_searcher.set_grid_search_center(sh2_);

		// Get the relevant spin coordinates from the pose
		single_dataset_tb.update_spin_coordinates(sh2_);

		// Run the gridsearch
		core::Real score;
		core::Real best_score(std::numeric_limits<core::Real>::max());
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords;
		while ( grid_searcher.valid_next_grid_point(metal_coords) ) {
			score = single_dataset_tb.solve_tensor_and_compute_score_by_svd(metal_coords) * weight;
			if ( score < best_score ) {
				best_score = score;
				best_metal_coords = metal_coords;
				grid_searcher.set_best_grid_point(best_metal_coords);
			}
		}

		// Best score and metal coordinates under given grid search parameters
		TS_ASSERT_DELTA(best_metal_coords.x(), -13.5235,1e-1);
		TS_ASSERT_DELTA(best_metal_coords.y(), -2.6035, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords.z(), -5.9535, 1e-1);
		score = single_dataset_tb.solve_tensor_and_compute_score_by_svd(best_metal_coords);

		// Show score and determined tensor
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tb.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << best_score << std::endl;
		TR.Debug << "PCS score = " << score << std::endl;
		PCSTensorOP tensor(single_dataset_tb.get_tensor());
		tensor->show_tensor_stats(TR.Debug);
		tensor->diagonalize_tensor();
		tensor->reorder_tensor();
		TS_ASSERT_DELTA(tensor->get_T_xx(),    -13.812, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_xy(),     -3.959, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_xz(),      3.872, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_yy(),      9.930, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_yz(),      7.090, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_xx(),   1.1438, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_yy(),  14.6657, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_zz(), -15.8095, 1e-1);
		TS_ASSERT_DELTA(tensor->get_ax(),      -23.714, 1e-1);
		TS_ASSERT_DELTA(tensor->get_rh(),      -13.522, 1e-1);
		TS_ASSERT_DELTA(tensor->get_alpha(),   12.9994, 1e-1);
		TS_ASSERT_DELTA(tensor->get_beta(),   105.2491, 1e-1);
		TS_ASSERT_DELTA(tensor->get_gamma(),   33.8767, 1e-1);
		tensor->show_tensor_stats(TR.Debug);
	}

	/// @brief test PCS tensor calculation for scaled data with grid search and SVD
	void test_calc_tensor_and_score_svd_on_scaled_data() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;
		using core::pose::named_atom_id_to_atom_id;

		core_init_with_additional_options("-nmr:pcs:normalize_data");

		TR << "Testing PCS tensor and score calculation by SVD" << std::endl;

		// create pcs dataset
		core::Real weight(1.0);
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, weight, "CONST", "SVD");

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearch grid_searcher(grid_atom1, grid_atom2, sh2_, 15, 1, 0, 25);

		// set grid search center coordinates
		grid_searcher.set_grid_search_center(sh2_);

		// Get the relevant spin coordinates from the pose
		single_dataset_tb.update_spin_coordinates(sh2_);

		// Run the gridsearch
		core::Real score;
		core::Real best_score(std::numeric_limits<core::Real>::max());
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords;
		while ( grid_searcher.valid_next_grid_point(metal_coords) ) {
			score = single_dataset_tb.solve_tensor_and_compute_score_by_svd(metal_coords) * weight;
			if ( score < best_score ) {
				best_score = score;
				best_metal_coords = metal_coords;
				grid_searcher.set_best_grid_point(best_metal_coords);
			}
		}

		// Best score and metal coordinates under given grid search parameters
		TS_ASSERT_DELTA(best_metal_coords.x(), -13.5235,1e-1);
		TS_ASSERT_DELTA(best_metal_coords.y(), -2.6035, 1e-1);
		TS_ASSERT_DELTA(best_metal_coords.z(), -5.9535, 1e-1);
		score = single_dataset_tb.solve_tensor_and_compute_score_by_svd(best_metal_coords);

		// Show score and determined tensor
		TR.Debug << "Calculated PCS score and tensor for (scaled) dataset " << single_dataset_tb.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << best_score << std::endl;
		TR.Debug << "PCS score = " << score << std::endl;
		PCSTensorOP tensor(single_dataset_tb.get_tensor());
		tensor->show_tensor_stats(TR.Debug);
		tensor->diagonalize_tensor();
		tensor->reorder_tensor();
		TS_ASSERT_DELTA(tensor->get_T_xx(),    -13.812, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_xy(),     -3.959, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_xz(),      3.872, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_yy(),      9.930, 1e-1);
		TS_ASSERT_DELTA(tensor->get_T_yz(),      7.090, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_xx(),   1.1438, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_yy(),  14.6657, 1e-1);
		TS_ASSERT_DELTA(tensor->get_Eig_zz(), -15.8095, 1e-1);
		TS_ASSERT_DELTA(tensor->get_ax(),      -23.714, 1e-1);
		TS_ASSERT_DELTA(tensor->get_rh(),      -13.522, 1e-1);
		TS_ASSERT_DELTA(tensor->get_alpha(),   12.9994, 1e-1);
		TS_ASSERT_DELTA(tensor->get_beta(),   105.2491, 1e-1);
		TS_ASSERT_DELTA(tensor->get_gamma(),   33.8767, 1e-1);
		tensor->show_tensor_stats(TR.Debug);
	}

	/// @brief test PCS tensor and score calculation with NLS
	void test_calc_tensor_and_score_nls() {
		using namespace core::scoring::nmr::pcs;

		TR << "Testing PCS tensor and score calculation by NLS" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// create pcs dataset
		core::Real weight(1.0);
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, weight, "CONST", "NLS");
		PCSSingleSet single_dataset_dy("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, weight, "CONST", "NLS");
		PCSSingleSet single_dataset_tm("core/scoring/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_, weight, "CONST", "NLS");
		PCSSingleSet single_dataset_er("core/scoring/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_, weight, "CONST", "NLS");

		// Get the relevant spin coordinates from the pose
		single_dataset_tb.update_spin_coordinates(sh2_);
		single_dataset_dy.update_spin_coordinates(sh2_);
		single_dataset_tm.update_spin_coordinates(sh2_);
		single_dataset_er.update_spin_coordinates(sh2_);

		// Set metal start coordinates and coordinate limits
		numeric::xyzVector< core::Real > start_metal_coords(-13.0, -3.0, -6.0);
		utility::fixedsizearray1<core::Real,6> metal_coords_bounds;
		metal_coords_bounds[1] = -23.0; metal_coords_bounds[2] = -13.0; metal_coords_bounds[3] = -16.0;
		metal_coords_bounds[4] =  -3.0; metal_coords_bounds[5] =   7.0; metal_coords_bounds[6] = 4.0;
		single_dataset_tb.set_metal_coord_bounds(metal_coords_bounds);
		single_dataset_dy.set_metal_coord_bounds(metal_coords_bounds);
		single_dataset_tm.set_metal_coord_bounds(metal_coords_bounds);
		single_dataset_er.set_metal_coord_bounds(metal_coords_bounds);

		// Run tensor and score calculation
		core::Real score_dataset_tb = single_dataset_tb.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_dy = single_dataset_dy.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_tm = single_dataset_tm.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_er = single_dataset_er.solve_tensor_and_compute_score_by_nls(start_metal_coords);

		// Show score and determined tensor Tb
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tb.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_tb << std::endl;
		PCSTensorOP tensor_tb_dataset(single_dataset_tb.get_tensor());
		tensor_tb_dataset->reorder_tensor();
		tensor_tb_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_ax(), -22.528, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_rh(), -10.775, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().x(), -12.993, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().y(),  -2.884, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().z(),  -5.876, 1e-1);

		// Show score and determined tensor Dy
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_dy.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_dy << std::endl;
		PCSTensorOP tensor_dy_dataset(single_dataset_dy.get_tensor());
		tensor_dy_dataset->reorder_tensor();
		tensor_dy_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_ax(), -17.600, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_rh(),  -8.230, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().x(), -12.996, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().y(),  -2.882, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().z(),  -5.877, 1e-1);

		// Show score and determined tensor Tm
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tm.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_tm << std::endl;
		PCSTensorOP tensor_tm_dataset(single_dataset_tm.get_tensor());
		tensor_tm_dataset->reorder_tensor();
		tensor_tm_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_ax(), 19.833, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_rh(),  3.193, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().x(), -12.994, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().y(),  -2.884, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().z(),  -5.876, 1e-1);

		// Show score and determined tensor Er
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_er.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_er << std::endl;
		PCSTensorOP tensor_er_dataset(single_dataset_er.get_tensor());
		tensor_er_dataset->reorder_tensor();
		tensor_er_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_er_dataset->get_ax(), 7.263, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_rh(), 2.274, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().x(), -12.991, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().y(),  -2.883, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().z(),  -5.875, 1e-1);

	}

	/// @brief test PCS tensor and score calculation with NLS
	///        and fixed predefined values for Xax and Xrh
	void test_calc_tensor_and_score_nls_fixed_XaxXrh() {
		using namespace core::scoring::nmr::pcs;

		TR << "Testing PCS tensor and score calculation by NLS with fixed Xax and Xrh" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// create pcs dataset
		core::Real weight(1.0);
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, weight, "CONST", "NLSAX");
		PCSSingleSet single_dataset_dy("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, weight, "CONST", "NLSAXRH");
		PCSSingleSet single_dataset_tm("core/scoring/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_, weight, "CONST", "NLSRH");
		PCSSingleSet single_dataset_er("core/scoring/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_, weight, "CONST", "NLSAXRH");

		// Get the relevant spin coordinates from the pose
		single_dataset_tb.update_spin_coordinates(sh2_);
		single_dataset_dy.update_spin_coordinates(sh2_);
		single_dataset_tm.update_spin_coordinates(sh2_);
		single_dataset_er.update_spin_coordinates(sh2_);

		// Set metal start coordinates and coordinate limits
		numeric::xyzVector< core::Real > start_metal_coords(-13.0, -3.0, -6.0);
		utility::fixedsizearray1<core::Real,6> metal_coords_bounds;
		metal_coords_bounds[1] = -23.0; metal_coords_bounds[2] = -13.0; metal_coords_bounds[3] = -16.0;
		metal_coords_bounds[4] =  -3.0; metal_coords_bounds[5] =   7.0; metal_coords_bounds[6] = 4.0;

		utility::vector1< core::Real > fixed_tensor_values(8);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360; // alpha, beta, gamma
		fixed_tensor_values[4] = fixed_tensor_values[5] = fixed_tensor_values[6] = 99.9; // xM, yM, zM
		fixed_tensor_values[7] = fixed_tensor_values[8] = 99.9; // Xax, Xrh

		single_dataset_tb.set_metal_coord_bounds(metal_coords_bounds);
		fixed_tensor_values[7] = -22.528;
		PCSTensorOP tb_tensor( new PCSTensor() );
		tb_tensor->set_tensor_in_pas(fixed_tensor_values);
		single_dataset_tb.set_tensor(tb_tensor);

		single_dataset_dy.set_metal_coord_bounds(metal_coords_bounds);
		fixed_tensor_values[7] = -17.569;
		fixed_tensor_values[8] = -8.225;
		PCSTensorOP dy_tensor( new PCSTensor() );
		dy_tensor->set_tensor_in_pas(fixed_tensor_values);
		single_dataset_dy.set_tensor(dy_tensor);

		single_dataset_tm.set_metal_coord_bounds(metal_coords_bounds);
		fixed_tensor_values[8] = 3.193;
		PCSTensorOP tm_tensor( new PCSTensor() );
		tm_tensor->set_tensor_in_pas(fixed_tensor_values);
		single_dataset_tm.set_tensor(tm_tensor);

		single_dataset_er.set_metal_coord_bounds(metal_coords_bounds);
		fixed_tensor_values[7] = 7.263;
		fixed_tensor_values[8] = 2.274;
		PCSTensorOP er_tensor( new PCSTensor() );
		er_tensor->set_tensor_in_pas(fixed_tensor_values);
		single_dataset_er.set_tensor(er_tensor);

		// Run tensor and score calculation
		core::Real score_dataset_tb = single_dataset_tb.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_dy = single_dataset_dy.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_tm = single_dataset_tm.solve_tensor_and_compute_score_by_nls(start_metal_coords);
		core::Real score_dataset_er = single_dataset_er.solve_tensor_and_compute_score_by_nls(start_metal_coords);

		// Show score and determined tensor Tb
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tb.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_tb << std::endl;
		PCSTensorOP tensor_tb_dataset(single_dataset_tb.get_tensor());
		tensor_tb_dataset->reorder_tensor();
		tensor_tb_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_ax(), -22.528, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_rh(), -10.775, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().x(), -12.993, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().y(),  -2.884, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().z(),  -5.876, 1e-1);

		// Show score and determined tensor Dy
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_dy.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_dy << std::endl;
		PCSTensorOP tensor_dy_dataset(single_dataset_dy.get_tensor());
		tensor_dy_dataset->reorder_tensor();
		tensor_dy_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_ax(), -17.569, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_rh(),  -8.230, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().x(), -12.996, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().y(),  -2.882, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().z(),  -5.877, 1e-1);

		// Show score and determined tensor Tm
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tm.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_tm << std::endl;
		PCSTensorOP tensor_tm_dataset(single_dataset_tm.get_tensor());
		tensor_tm_dataset->reorder_tensor();
		tensor_tm_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_ax(), 19.8807, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_rh(),   3.193, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().x(), -13.0055, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().y(),  -2.9358, 1e-1);
		TS_ASSERT_DELTA(tensor_tm_dataset->get_metal_center().z(),  -5.8207, 1e-1);

		// Show score and determined tensor Er
		TR.Debug << "Best calculated PCS score and tensor for dataset " << single_dataset_er.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_er << std::endl;
		PCSTensorOP tensor_er_dataset(single_dataset_er.get_tensor());
		tensor_er_dataset->reorder_tensor();
		tensor_er_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_er_dataset->get_ax(), 7.263, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_rh(), 2.274, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().x(), -12.990, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().y(),  -2.882, 1e-1);
		TS_ASSERT_DELTA(tensor_er_dataset->get_metal_center().z(),  -5.875, 1e-1);

	}

	/// @brief test PCS tensor and score calculation with NLS and SVD
	///        and with different weighting schemes
	void test_calc_tensor_and_score_nls_svd_different_weightings() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		TR << "Testing PCS tensor and score calculation using different weighting methods" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// create pcs dataset
		core::Real weight(1.0);
		PCSSingleSet single_dataset_tb("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, weight, "SIGMA", "SVD");
		PCSSingleSet single_dataset_dy("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, weight, "OBSIG", "NLS");

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearch grid_searcher(grid_atom1, grid_atom2, sh2_, 15, 1, 0, 25);

		// set grid search center coordinates
		grid_searcher.set_grid_search_center(sh2_);

		// Get the relevant spin coordinates from the pose
		single_dataset_tb.update_spin_coordinates(sh2_);
		single_dataset_dy.update_spin_coordinates(sh2_);

		// Set metal start coordinates and coordinate limits for NLS fitting
		numeric::xyzVector< core::Real > start_metal_coords(-13.0, -3.0, -6.0);
		utility::fixedsizearray1<core::Real,6> metal_coords_bounds;
		metal_coords_bounds[1] = -23.0; metal_coords_bounds[2] = -13.0; metal_coords_bounds[3] = -16.0;
		metal_coords_bounds[4] =  -3.0; metal_coords_bounds[5] =   7.0; metal_coords_bounds[6] = 4.0;
		single_dataset_dy.set_metal_coord_bounds(metal_coords_bounds);

		// Run the gridsearch and calculate tensor by SVD
		core::Real score_dataset_tb;
		core::Real best_score_dataset_tb(std::numeric_limits<core::Real>::max());
		numeric::xyzVector<core::Real> metal_coords;
		numeric::xyzVector<core::Real> best_metal_coords;
		while ( grid_searcher.valid_next_grid_point(metal_coords) ) {
			score_dataset_tb = single_dataset_tb.solve_tensor_and_compute_score_by_svd(metal_coords) * weight;
			if ( score_dataset_tb < best_score_dataset_tb ) {
				best_score_dataset_tb = score_dataset_tb;
				best_metal_coords = metal_coords;
				grid_searcher.set_best_grid_point(best_metal_coords);
			}
		}
		score_dataset_tb = single_dataset_tb.solve_tensor_and_compute_score_by_svd(best_metal_coords);

		// Run tensor and score calculation by NLS fitting
		core::Real score_dataset_dy = single_dataset_dy.solve_tensor_and_compute_score_by_nls(start_metal_coords);

		// Show score and tensor as determined by SVD for Tb dataset
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_tb.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << best_score_dataset_tb << std::endl;
		TR.Debug << "PCS score = " << score_dataset_tb << std::endl;
		PCSTensorOP tensor_tb_dataset(single_dataset_tb.get_tensor());
		tensor_tb_dataset->show_tensor_stats(TR.Debug);
		tensor_tb_dataset->diagonalize_tensor();
		tensor_tb_dataset->reorder_tensor();
		tensor_tb_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_T_xx(),    -13.812, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_T_xy(),     -3.959, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_T_xz(),      3.872, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_T_yy(),      9.930, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_T_yz(),      7.090, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_Eig_xx(),   1.1438, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_Eig_yy(),  14.6657, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_Eig_zz(), -15.8095, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_ax(),      -23.714, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_rh(),      -13.522, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_alpha(),   12.9994, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_beta(),   105.2491, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_gamma(),   33.8767, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().x(), -13.5235, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().y(),  -2.6035, 1e-1);
		TS_ASSERT_DELTA(tensor_tb_dataset->get_metal_center().z(),  -5.9535, 1e-1);

		// Show score and tensor as determined by NLS fitting for Dy dataset
		TR.Debug << "Calculated PCS score and tensor for dataset " << single_dataset_dy.get_dataset_name() << std::endl;
		TR.Debug << "PCS score = " << score_dataset_dy << std::endl;
		PCSTensorOP tensor_dy_dataset(single_dataset_dy.get_tensor());
		tensor_dy_dataset->reorder_tensor();
		tensor_dy_dataset->show_tensor_stats(TR.Debug);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_ax(), -17.603, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_rh(),  -8.230, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().x(), -12.995, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().y(),  -2.883, 1e-1);
		TS_ASSERT_DELTA(tensor_dy_dataset->get_metal_center().z(),  -5.875, 1e-1);
	}
};
