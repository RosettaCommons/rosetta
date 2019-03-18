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
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
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

static basic::Tracer TR("core.scoring.nmr.pcs.PCSMultiSet.cxxtest");

class PCSMultiSetTests : public CxxTest::TestSuite {

private:
	core::pose::Pose sh2_;
	core::pose::Pose cam_;

public:

	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(sh2_, "core/scoring/nmr/pcs/1x0n_model1_renum.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		sh2_.clear();
		cam_.clear();
	}

	/// @brief test creation of PCSMultiSet
	void test_PCSMultiSet_instantiation_and_getters() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		TR << "Test PCSMultiSet instantiation with NMRGridSearch and NMRSpinlabel" << std::endl;

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearchOP gridsearch( new NMRGridSearch(grid_atom1, grid_atom2, sh2_, 15, 1, 0, 25) );

		// create spinlabel object
		NMRSpinlabelOP spinlabel( new NMRSpinlabel("fa_standard", "R1A") );

		// Create PCSMultiSet with gridsearch
		utility::vector1<PCSSingleSetOP> pcs_data_sh2;
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, 1.0, "CONST", "SVD")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, 1.0, "CONST", "NLS")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_, 1.0, "CONST", "SVD")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_, 1.0, "CONST", "NLS")));
		PCSMultiSetOP pcs_sh2( new PCSMultiSet(pcs_data_sh2, 95, gridsearch) );

		// Create PCSMultiSet with spinlabel
		utility::vector1<PCSSingleSetOP> pcs_data_cam;
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tb.txt", "Tb", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_dy.txt", "Dy", cam_, 1.0, "CONST", "NLS")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tm.txt", "Tm", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_yb.txt", "Yb", cam_, 1.0, "CONST", "NLS")));
		PCSMultiSetOP pcs_cam( new PCSMultiSet(pcs_data_cam, 13, spinlabel) );

		TS_ASSERT_EQUALS(pcs_sh2->get_number_metal_ions(), 4);
		TS_ASSERT_EQUALS(pcs_cam->get_number_metal_ions(), 4);
		TS_ASSERT_EQUALS(pcs_sh2->get_pcs_singleset_vec().size(), 4);
		TS_ASSERT_EQUALS(pcs_cam->get_pcs_singleset_vec().size(), 4);
		TS_ASSERT_EQUALS(pcs_sh2->get_tag_residue_number(), 95);
		TS_ASSERT_EQUALS(pcs_cam->get_tag_residue_number(), 13);
		TS_ASSERT(pcs_sh2->get_gridsearch_iterator());
		TS_ASSERT_EQUALS(pcs_sh2->get_gridsearch_iterator()->get_grid_atom1().rsd(), 95);
		TS_ASSERT_DELTA(pcs_sh2->get_gridsearch_iterator()->get_grid_search_center().x(), -10.5235832, 1.0e-2);
		TS_ASSERT_DELTA(pcs_sh2->get_gridsearch_iterator()->get_grid_search_center().y(), -13.6039896, 1.0e-2);
		TS_ASSERT_DELTA(pcs_sh2->get_gridsearch_iterator()->get_grid_search_center().z(),  -3.9542288, 1.0e-2);
		TS_ASSERT(pcs_cam->get_spinlabel());
		TS_ASSERT_EQUALS(pcs_cam->get_spinlabel()->get_code(), "R1A");
		TS_ASSERT_EQUALS(pcs_cam->get_spinlabel()->get_radical_atom(), "O1");

	}

	/// @brief test score calculation of PCSMultiSet with grid search
	void test_PCSMultiSet_score_calculation_with_gridsearch() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		TR << "Test PCSMultiSet score calculation with NMRGridSearch" << std::endl;
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		// create gridsearch object
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("CA", 95), sh2_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("CB", 95), sh2_) );
		NMRGridSearchOP gridsearch( new NMRGridSearch(grid_atom1, grid_atom2, sh2_, 15, 1, 0, 25) );

		// Create PCSMultiSet with gridsearch
		utility::vector1<PCSSingleSetOP> pcs_data_sh2;
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_tb.txt", "Tb", sh2_, 1.0, "CONST", "SVD")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_dy.txt", "Dy", sh2_, 1.0, "CONST", "NLS")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_tm.txt", "Tm", sh2_, 1.0, "CONST", "SVD")));
		pcs_data_sh2.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_er.txt", "Er", sh2_, 1.0, "CONST", "NLS")));
		PCSMultiSetOP pcs_sh2( new PCSMultiSet(pcs_data_sh2, 95, gridsearch) );

		utility::vector1<PCSTensorCOP> tensors;
		core::Real total_score = pcs_sh2->compute_score(sh2_, tensors);
		TR.Debug << "Total PCS score is " << total_score << std::endl;

		PCSTensorOP tensor_tb( new PCSTensor( *tensors[1] ) );
		tensor_tb->diagonalize_tensor();
		tensor_tb->reorder_tensor();
		TS_ASSERT_DELTA(tensor_tb->get_T_xx(),    -13.812, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_xy(),     -3.959, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_xz(),      3.872, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_yy(),      9.930, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_yz(),      7.090, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_ax(),      -23.714, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_rh(),      -13.522, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().x(), -13.524, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().y(),  -2.604, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().z(),  -5.954, 1.e-1);

		PCSTensorOP tensor_dy( new PCSTensor( *tensors[2] ) );
		tensor_dy->reorder_tensor();
		TS_ASSERT_DELTA(tensor_dy->get_ax(),    -17.600, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_rh(),     -8.230, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().x(), -12.996, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().y(),  -2.882, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().z(),  -5.877, 1.e-1);

		PCSTensorOP tensor_tm( new PCSTensor( *tensors[3] ) );
		tensor_tm->diagonalize_tensor();
		tensor_tm->reorder_tensor();
		TS_ASSERT_DELTA(tensor_tm->get_T_xx(),    12.179, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_xy(),     1.735, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_xz(),    -5.180, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_yy(),    -8.968, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_yz(),    -0.330, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_ax(),     20.8418, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_rh(),      4.3271, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().x(), -13.524, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().y(),  -2.604, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().z(),  -5.954, 1.e-1);

		PCSTensorOP tensor_er( new PCSTensor( *tensors[4] ) );
		tensor_er->reorder_tensor();
		TS_ASSERT_DELTA(tensor_er->get_ax(),    7.263, 1.e-1);
		TS_ASSERT_DELTA(tensor_er->get_rh(),    2.274, 1.e-1);
		TS_ASSERT_DELTA(tensor_er->get_metal_center().x(), -12.991, 1.e-1);
		TS_ASSERT_DELTA(tensor_er->get_metal_center().y(),  -2.883, 1.e-1);
		TS_ASSERT_DELTA(tensor_er->get_metal_center().z(),  -5.875, 1.e-1);

	}

	/// @brief test score calculation of PCSMultiSet with spinlabel
	void test_PCSMultiSet_score_calculation_with_spinlabel() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		core_init_with_additional_options("-nmr:spinlabel:highres_conformer_filter_type DISTANCE");
		// Set fixed RG seed for NLS fitting
		initialize_rng();

		TR << "Test PCSMultiSet score calculation with NMRSpinlabel" << std::endl;

		// create spinlabel object
		NMRSpinlabelOP spinlabel( new NMRSpinlabel("fa_standard", "R1A") );

		// Create PCSMultiSet with spinlabel
		utility::vector1<PCSSingleSetOP> pcs_data_cam;
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tb.txt", "Tb", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_dy.txt", "Dy", cam_, 1.0, "CONST", "NLS")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tm.txt", "Tm", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_yb.txt", "Yb", cam_, 1.0, "CONST", "NLS")));
		PCSMultiSetOP pcs_cam( new PCSMultiSet(pcs_data_cam, 13, spinlabel) );

		utility::vector1<PCSTensorCOP> tensors;
		core::Real total_score = pcs_cam->compute_score(cam_, tensors);
		TR.Debug << "Total PCS score is " << total_score << std::endl;

		PCSTensorOP tensor_tb( new PCSTensor( *tensors[1] ) );
		tensor_tb->diagonalize_tensor();
		tensor_tb->reorder_tensor();
		TS_ASSERT_DELTA(tensor_tb->get_T_xx(),    3.4702, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_xy(),   -1.8080, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_xz(),   -2.0265, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_yy(),    6.7539, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_T_yz(),   -2.2196, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_ax(),    -16.2946, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_rh(),     -4.4245, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().x(), 85.9314, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().y(), 32.8311, 1.e-1);
		TS_ASSERT_DELTA(tensor_tb->get_metal_center().z(),  9.3183, 1.e-1);

		PCSTensorOP tensor_dy( new PCSTensor( *tensors[2] ) );
		tensor_dy->reorder_tensor();
		TS_ASSERT_DELTA(tensor_dy->get_ax(),    -25.0, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_rh(),     -8.0, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().x(), 86.679, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().y(), 33.309, 1.e-1);
		TS_ASSERT_DELTA(tensor_dy->get_metal_center().z(), 10.941, 1.e-1);

		PCSTensorOP tensor_tm( new PCSTensor( *tensors[3] ) );
		tensor_tm->diagonalize_tensor();
		tensor_tm->reorder_tensor();
		TS_ASSERT_DELTA(tensor_tm->get_T_xx(),   -2.5307, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_xy(),    1.5964, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_xz(),    1.4392, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_yy(),   -4.2935, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_T_yz(),    1.3550, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_ax(),     10.8894, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_rh(),      3.2460, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().x(), 85.9314, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().y(), 32.8311, 1.e-1);
		TS_ASSERT_DELTA(tensor_tm->get_metal_center().z(),  9.3183, 1.e-1);

		PCSTensorOP tensor_yb( new PCSTensor( *tensors[4] ) );
		tensor_yb->reorder_tensor();
		TS_ASSERT_DELTA(tensor_yb->get_ax(),      8.0, 1.e-1);
		TS_ASSERT_DELTA(tensor_yb->get_rh(),      2.0, 1.e-1);
		TS_ASSERT_DELTA(tensor_yb->get_metal_center().x(), 86.679, 1.e-1);
		TS_ASSERT_DELTA(tensor_yb->get_metal_center().y(), 33.309, 1.e-1);
		TS_ASSERT_DELTA(tensor_yb->get_metal_center().z(), 10.941, 1.e-1);

	}

	/// @brief test score calculation of PCSMultiSet with fixed input tensor
	void test_PCSMultiSet_score_calculation_from_fixed_tensor() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pcs;

		TR << "Test PCSMultiSet score calculation with fixed input tensor" << std::endl;

		// Create PCSMultiSet with default spinlabel
		utility::vector1<PCSSingleSetOP> pcs_data_cam;
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tb.txt", "Tb", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_dy.txt", "Dy", cam_, 1.0, "CONST", "NLS")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_tm.txt", "Tm", cam_, 1.0, "CONST", "SVD")));
		pcs_data_cam.emplace_back(PCSSingleSetOP( new PCSSingleSet("core/scoring/nmr/pcs/sim_pcs_val_mtsl_13_yb.txt", "Yb", cam_, 1.0, "CONST", "NLS")));
		PCSMultiSetOP pcs_cam( new PCSMultiSet(pcs_data_cam, 13) );

		// Set fixed tensor
		utility::vector1< core::Real > fixed_tensor_values(8);
		fixed_tensor_values[1] = fixed_tensor_values[2] = fixed_tensor_values[3] = 360; // alpha, beta, gamma
		fixed_tensor_values[4] = fixed_tensor_values[5] = fixed_tensor_values[6] = 99.9; // xM, yM, zM
		fixed_tensor_values[7] = fixed_tensor_values[8] = 99.9; // Xax, Xrh

		PCSTensorOP tensor_tb( new PCSTensor() );
		fixed_tensor_values[1] = 10.0;   fixed_tensor_values[2] = 10.0;   fixed_tensor_values[3] = 10.0;
		fixed_tensor_values[4] = 86.679; fixed_tensor_values[5] = 33.309; fixed_tensor_values[6] = 10.941;
		fixed_tensor_values[7] = -20.0;  fixed_tensor_values[8] = -4.0;
		tensor_tb->set_tensor_in_pas(fixed_tensor_values);
		pcs_cam->get_pcs_singleset_vec()[1]->set_tensor(tensor_tb);

		PCSTensorOP tensor_dy( new PCSTensor() );
		fixed_tensor_values[1] = 12.0;   fixed_tensor_values[2] =  8.0;   fixed_tensor_values[3] = 12.0;
		fixed_tensor_values[4] = 86.679; fixed_tensor_values[5] = 33.309; fixed_tensor_values[6] = 10.941;
		fixed_tensor_values[7] = -25.0;  fixed_tensor_values[8] = -8.0;
		tensor_dy->set_tensor_in_pas(fixed_tensor_values);
		pcs_cam->get_pcs_singleset_vec()[2]->set_tensor(tensor_dy);

		PCSTensorOP tensor_tm( new PCSTensor() );
		fixed_tensor_values[1] = 15.0;   fixed_tensor_values[2] = 10.0;   fixed_tensor_values[3] = 12.0;
		fixed_tensor_values[4] = 86.679; fixed_tensor_values[5] = 33.309; fixed_tensor_values[6] = 10.941;
		fixed_tensor_values[7] = 15.0;   fixed_tensor_values[8] = 4.0;
		tensor_tm->set_tensor_in_pas(fixed_tensor_values);
		pcs_cam->get_pcs_singleset_vec()[3]->set_tensor(tensor_tm);

		PCSTensorOP tensor_yb( new PCSTensor() );
		fixed_tensor_values[1] =  5.0;   fixed_tensor_values[2] =  5.0;   fixed_tensor_values[3] = 10.0;
		fixed_tensor_values[4] = 86.679; fixed_tensor_values[5] = 33.309; fixed_tensor_values[6] = 10.941;
		fixed_tensor_values[7] =  8.0;   fixed_tensor_values[8] = 2.0;
		tensor_yb->set_tensor_in_pas(fixed_tensor_values);
		pcs_cam->get_pcs_singleset_vec()[4]->set_tensor(tensor_yb);

		// Set tensor for this multiset fixed
		pcs_cam->fix_tensors();

		utility::vector1<PCSTensorCOP> tensors;
		core::Real total_score = pcs_cam->compute_score(cam_, tensors);
		TR.Debug << "Total PCS score is " << total_score << std::endl;

		// Score should be close to zero
		TS_ASSERT_DELTA(total_score, 0.0, 1.0e-3);

		// Predicted PCSs must be close to input PCSs
		// Tb3+
		const core::Real tol(1.0e-2);
		utility::vector1<PCSSingle> const & pcs_tb = pcs_cam->get_pcs_singleset_vec()[1]->get_single_pcs_vec();
		TS_ASSERT_EQUALS(pcs_tb.size(),124);
		TS_ASSERT_DELTA(pcs_tb[5].get_pcs_calc(), -0.866, tol);
		TS_ASSERT_DELTA(pcs_tb[18].get_pcs_calc(), 0.558, tol);
		TS_ASSERT_DELTA(pcs_tb[42].get_pcs_calc(), 0.681, tol);
		TS_ASSERT_DELTA(pcs_tb[51].get_pcs_calc(),-0.529, tol);
		TS_ASSERT_DELTA(pcs_tb[93].get_pcs_calc(), 1.979, tol);

		// Dy3+
		utility::vector1<PCSSingle> const & pcs_dy = pcs_cam->get_pcs_singleset_vec()[2]->get_single_pcs_vec();
		TS_ASSERT_EQUALS(pcs_dy.size(),123);
		TS_ASSERT_DELTA(pcs_dy[5].get_pcs_calc(),  1.119, tol);
		TS_ASSERT_DELTA(pcs_dy[18].get_pcs_calc(), 0.493, tol);
		TS_ASSERT_DELTA(pcs_dy[42].get_pcs_calc(), 0.549, tol);
		TS_ASSERT_DELTA(pcs_dy[51].get_pcs_calc(),-0.482, tol);
		TS_ASSERT_DELTA(pcs_dy[93].get_pcs_calc(), 0.889, tol);

		// Tm3+
		utility::vector1<PCSSingle> const & pcs_tm = pcs_cam->get_pcs_singleset_vec()[3]->get_single_pcs_vec();
		TS_ASSERT_EQUALS(pcs_tm.size(),127);
		TS_ASSERT_DELTA(pcs_tm[5].get_pcs_calc(),  1.860, tol);
		TS_ASSERT_DELTA(pcs_tm[18].get_pcs_calc(),-0.655, tol);
		TS_ASSERT_DELTA(pcs_tm[42].get_pcs_calc(),-0.237, tol);
		TS_ASSERT_DELTA(pcs_tm[51].get_pcs_calc(), 0.298, tol);
		TS_ASSERT_DELTA(pcs_tm[93].get_pcs_calc(),-0.344, tol);

		// Yb3+
		utility::vector1<PCSSingle> const & pcs_yb = pcs_cam->get_pcs_singleset_vec()[4]->get_single_pcs_vec();
		TS_ASSERT_EQUALS(pcs_yb.size(),136);
		TS_ASSERT_DELTA(pcs_yb[5].get_pcs_calc(),  1.124, tol);
		TS_ASSERT_DELTA(pcs_yb[18].get_pcs_calc(),-1.038, tol);
		TS_ASSERT_DELTA(pcs_yb[42].get_pcs_calc(),-0.091, tol);
		TS_ASSERT_DELTA(pcs_yb[51].get_pcs_calc(),-0.106, tol);
		TS_ASSERT_DELTA(pcs_yb[93].get_pcs_calc(),-0.043, tol);
	}

};
