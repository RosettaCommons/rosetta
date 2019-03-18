// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREMultiSet.cxxtest.hh
/// @brief   unit test for class PREMultiSet that stores and handles data of multiple PRE experiments collected for the same spinlabel
/// @details Last modified: 10/16/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/util.hh>
#include <core/io/nmr/ParaIon.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

static basic::Tracer TR("core.scoring.nmr.pre.PREMultiSet.cxxtest");

class PREMultiSetTests : public CxxTest::TestSuite {

public:
	using WeightCoordVector = core::scoring::nmr::pre::PREMultiSet::WeightCoordVector;

private:
	core::pose::Pose cam_;
	core::pose::Pose mid1_;
	WeightCoordVector radical_atom_wts_coords_;

public:

	/// @brief Setup Test
	void setUp() {

		// Initialize core & options system
		core_init();

		// Load pose from pdb
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(mid1_, "core/scoring/nmr/pre/3v1e_AB_Mn.pdb", core::import_pose::PDB_file);

		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(81.299, 33.311,  7.942)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(92.649, 31.000,  8.182)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(91.389, 32.009, 10.363)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.609, 34.451,  9.000)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.679, 33.309, 10.941)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(87.027, 35.888,  7.247)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(87.214, 35.445,  4.572)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.878, 34.328,  7.101)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.093, 32.216, 10.080)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.145, 32.886, 10.977)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(83.631, 34.837,  7.928)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.025, 35.951,  7.301)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(84.344, 30.100, 10.762)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.231, 32.869,  9.176)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.496, 33.889,  7.204)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.633, 30.224,  9.316)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.796, 34.025,  8.689)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(84.720, 29.332, 11.161)));
		radical_atom_wts_coords_.push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.801, 33.438,  7.847)));

	}

	void tearDown() {
		cam_.clear();
		mid1_.clear();
		radical_atom_wts_coords_.clear();
	}

	/// @brief test construction of PREMultiSet from data files and PRESingleSet vector
	///        and access of its data members
	void test_PREMultiSet_instantiation_and_getters_and_setters() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pre;

		TR << "Testing PREMultiSet creation" << std::endl;

		// Create first PREMultiSet with default spinlabel
		PRESingleSetOP pre_data_c13z_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_1h.txt", cam_, 1.0, "R2", "SIGMA"));
		PRESingleSetOP pre_data_c13z_r2_15n( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_15n.txt", cam_, 10.0, "R2", "CONST"));
		utility::vector1<PRESingleSetOP> singleset_vec;
		singleset_vec.push_back(pre_data_c13z_r2_1h);
		singleset_vec.push_back(pre_data_c13z_r2_15n);
		PREMultiSet pre_data_c13z(singleset_vec, cam_, 13, "Nitroxide");

		// Get number of PREs
		TS_ASSERT_EQUALS(pre_data_c13z_r2_1h->get_number_pre(), 132);
		TS_ASSERT_EQUALS(pre_data_c13z_r2_15n->get_number_pre(), 133);
		TS_ASSERT_EQUALS(pre_data_c13z.get_total_number_pre(), 265);

		// Get number of experiments
		TS_ASSERT_EQUALS(pre_data_c13z.get_number_experiments(), 2);

		// Access PRE data
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_pre_single_vec()[1].get_pre_exp(), 20.249, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z_r2_15n->get_pre_single_vec()[12].get_pre_exp(), 2.9015, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_values()[1], 20.249, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_values()[133], 0.2707, 1.0e-3);

		// Access PRE errors and weights
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_pre_single_vec()[1].get_pre_err(), 2.025, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z_r2_15n->get_pre_single_vec()[12].get_pre_err(), 0.2901, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_single_weights()[1], 0.2438, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_single_weights()[133], 1.000, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_weight(), 1.0, 1.0e-6);
		TS_ASSERT_DELTA(pre_data_c13z_r2_15n->get_weight(), 10.0, 1.0e-6);
		TS_ASSERT_DELTA(pre_data_c13z.get_weight(), 1.0, 1.0e-6);

		// Get rate type
		TS_ASSERT_EQUALS(pre_data_c13z_r2_1h->get_pre_rate_type(), R2_PARA);
		TS_ASSERT_EQUALS(pre_data_c13z_r2_15n->get_pre_rate_type(), R2_PARA);

		// Get single value weighting scheme
		TS_ASSERT_EQUALS(pre_data_c13z_r2_1h->get_single_pre_weighting_scheme(), SIGMA);
		TS_ASSERT_EQUALS(pre_data_c13z_r2_15n->get_single_pre_weighting_scheme(), CONST);

		// Set experimental parameters and access them
		pre_data_c13z_r2_1h->set_field_strength(600.0);
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_field_strength(), 600.0, 1.0e-6);

		pre_data_c13z.set_protein_mass(16.1);
		TS_ASSERT_DELTA(pre_data_c13z.get_protein_mass(), 16.1, 1.0e-6);

		pre_data_c13z.set_temperature(298.0);
		TS_ASSERT_DELTA(pre_data_c13z.get_temperature(), 298.0, 1.0e-6);

		pre_data_c13z.set_tau_c(1.0e-9);
		TS_ASSERT_DELTA(pre_data_c13z.get_tau_c(), 1.0e-9, 1.0e-10);

		// Ion type
		TS_ASSERT_EQUALS(pre_data_c13z.get_ion_type()->get_ion_label(), "Nitroxide");
		TS_ASSERT_DELTA(pre_data_c13z.get_ion_type()->get_tau_e()*1.0e-12, 100.0e-9, 1.0e-12);

		// Set and get spinlabel
		NMRSpinlabelOP spinlabel( new NMRSpinlabel("fa_standard", "R1A") );
		pre_data_c13z.set_spinlabel(spinlabel);
		TS_ASSERT_EQUALS(pre_data_c13z.get_spinlabel()->get_code(), "R1A");
		TS_ASSERT_EQUALS(pre_data_c13z.get_spinlabel_site_rsd(), 13);

		// Create another PREMultiSet for a metal protein and define grid search around metal site.
		// Grid search atom1 is the MN atom and grid search atom2 is the NE2 atom of the coordinating His15
		// Set the distance between atom1 and the center of the grid search to 0.0. In this way, the center
		// of the grid search box is in the MN atom.
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("MN",  44), mid1_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("NE2", 15), mid1_) );
		NMRGridSearchOP gridsearch( new NMRGridSearch(grid_atom1, grid_atom2, mid1_, 0, 0.5, 0, 1) );
		utility::vector1<PRESingleSetOP> singleset_vec_;
		PRESingleSetOP pres_mid1_mn_r1_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r1.txt", mid1_, 10.0, "R1", "CONST"));
		PRESingleSetOP pres_mid1_mn_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r2.txt", mid1_, 1.0, "R2", "CONST"));
		singleset_vec_.push_back(pres_mid1_mn_r1_1h);
		singleset_vec_.push_back(pres_mid1_mn_r2_1h);
		PREMultiSet pre_data_mid1_mn(singleset_vec_, mid1_, 44, "Mn2+", gridsearch);

		// Get number of PREs
		TS_ASSERT_EQUALS(pres_mid1_mn_r1_1h->get_number_pre(), 82);
		TS_ASSERT_EQUALS(pres_mid1_mn_r2_1h->get_number_pre(), 58);
		TS_ASSERT_EQUALS(pre_data_mid1_mn.get_total_number_pre(), 140);

		// Get number of experiments
		TS_ASSERT_EQUALS(pre_data_mid1_mn.get_number_experiments(), 2);

		// Access PRE data
		TS_ASSERT_DELTA(pres_mid1_mn_r1_1h->get_pre_single_vec()[10].get_pre_exp(), 166.990, 1.0e-3);
		TS_ASSERT_DELTA(pres_mid1_mn_r2_1h->get_pre_single_vec()[10].get_pre_exp(),  86.710, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_pre_values()[1], 0.480, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_pre_values()[83], 56.840, 1.0e-3);

		// Access PRE errors and weights
		TS_ASSERT_DELTA(pres_mid1_mn_r1_1h->get_pre_single_vec()[10].get_pre_err(), 10.00, 1.0e-3);
		TS_ASSERT_DELTA(pres_mid1_mn_r2_1h->get_pre_single_vec()[10].get_pre_err(), 25.00, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_pre_single_weights()[1], 1.00, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_pre_single_weights()[83], 1.00, 1.0e-3);
		TS_ASSERT_DELTA(pres_mid1_mn_r1_1h->get_weight(), 10.0, 1.0e-6);
		TS_ASSERT_DELTA(pres_mid1_mn_r2_1h->get_weight(), 1.0, 1.0e-6);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_weight(), 1.0, 1.0e-6);

		// Get rate type
		TS_ASSERT_EQUALS(pres_mid1_mn_r1_1h->get_pre_rate_type(), R1_PARA);
		TS_ASSERT_EQUALS(pres_mid1_mn_r2_1h->get_pre_rate_type(), R2_PARA);

		// Get single value weighting scheme
		TS_ASSERT_EQUALS(pres_mid1_mn_r1_1h->get_single_pre_weighting_scheme(), CONST);
		TS_ASSERT_EQUALS(pres_mid1_mn_r2_1h->get_single_pre_weighting_scheme(), CONST);

		// Set experimental parameters and access them
		pres_mid1_mn_r1_1h->set_field_strength(600.0);
		TS_ASSERT_DELTA(pres_mid1_mn_r1_1h->get_field_strength(), 600.0, 1.0e-6);

		pre_data_mid1_mn.set_protein_mass(10.2);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_protein_mass(), 10.2, 1.0e-6);

		pre_data_mid1_mn.set_temperature(298.0);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_temperature(), 298.0, 1.0e-6);

		pre_data_mid1_mn.set_tau_c(1.0e-9);
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_tau_c(), 1.0e-9, 1.0e-10);

		// Ion type
		TS_ASSERT_EQUALS(pre_data_mid1_mn.get_ion_type()->get_ion_label(), "Mn2+");
		TS_ASSERT_DELTA(pre_data_mid1_mn.get_ion_type()->get_tau_e()*1.0e-12, 5.0e-9, 1.0e-12);

	}

	/// @brief test construction of PREMultiSet from data files and PRESingleSet vector
	///        Apply normalization of input data by experiment StdDev
	void test_PREMultiSet_instantiation_and_data_scaling() {
		using namespace core::scoring::nmr::pre;
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		PRESingleSetOP pre_data_c13z_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_1h.txt", cam_, 1.0, "R2", "SIGMA"));
		PRESingleSetOP pre_data_c13z_r2_15n( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_15n.txt", cam_, 1.0, "R2", "OBSIG"));
		utility::vector1<PRESingleSetOP> singleset_vec;
		singleset_vec.push_back(pre_data_c13z_r2_1h);
		singleset_vec.push_back(pre_data_c13z_r2_15n);
		PREMultiSet pre_data_c13z(singleset_vec, cam_, 13, "Nitroxide");

		// Access single PRE values
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_pre_single_vec()[1].get_pre_exp(), 0.31915887, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z_r2_15n->get_pre_single_vec()[1].get_pre_exp(), 0.32820647, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_values()[1], 0.31915887, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_values()[133], 0.32820647, 1.0e-3);

		// Access single PRE errors and weights
		TS_ASSERT_DELTA(pre_data_c13z_r2_1h->get_pre_single_vec()[1].get_pre_err(), 0.03191746, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z_r2_15n->get_pre_single_vec()[1].get_pre_err(), 0.03285702, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_single_weights()[1], 981.61968569, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_c13z.get_pre_single_weights()[133], 59.9466336, 1.0e-3);

		TS_ASSERT(pre_data_c13z_r2_1h->normalized_data());
		TS_ASSERT(pre_data_c13z_r2_15n->normalized_data());

		utility::vector1<PRESingleSetOP> singleset_vec_;
		PRESingleSetOP pres_mid1_cu_r1_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_cu2_1h_r1.txt", mid1_, 10.0, "R1", "SIGMA"));
		PRESingleSetOP pres_mid1_cu_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_cu2_1h_r2.txt", mid1_, 1.0, "R2", "OBSIG"));
		singleset_vec_.push_back(pres_mid1_cu_r1_1h);
		singleset_vec_.push_back(pres_mid1_cu_r2_1h);
		PREMultiSet pre_data_mid1_cu(singleset_vec_, mid1_, 44, "Cu2+");

		TS_ASSERT_DELTA(pres_mid1_cu_r1_1h->get_pre_single_vec()[10].get_pre_exp(), 2.06018958, 1.0e-3);
		TS_ASSERT_DELTA(pres_mid1_cu_r2_1h->get_pre_single_vec()[10].get_pre_exp(), 5.03218813, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_cu.get_pre_values()[10], 2.06018958, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_cu.get_pre_values()[92], 5.03218813, 1.0e-3);

		TS_ASSERT_DELTA(pres_mid1_cu_r1_1h->get_pre_single_vec()[10].get_pre_err(), 0.27359755, 1.0e-3);
		TS_ASSERT_DELTA(pres_mid1_cu_r2_1h->get_pre_single_vec()[10].get_pre_err(), 0.12348019, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_cu.get_pre_single_weights()[10], 13.35905091, 1.0e-3);
		TS_ASSERT_DELTA(pre_data_mid1_cu.get_pre_single_weights()[92], 65.58513576, 1.0e-3);

		TS_ASSERT(pres_mid1_cu_r1_1h->normalized_data());
		TS_ASSERT(pres_mid1_cu_r2_1h->normalized_data());

	}

	/// @brief test fitting of spinlabel correlation time and R2 PRE score calculation
	void test_pre_r2_score_calculation_with_spinlabel() {
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		TR << "Testing PREMultiSet score calculation with spinlabel" << std::endl;

		PRESingleSetOP pre_data_c13z_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_1h.txt", cam_, 1.0, "R2", "CONST"));
		PRESingleSetOP pre_data_c13z_r2_15n( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_15n.txt", cam_, 1.0, "R2", "CONST"));
		utility::vector1<PRESingleSetOP> singleset_vec_c13z;
		singleset_vec_c13z.push_back(pre_data_c13z_r2_1h);
		singleset_vec_c13z.push_back(pre_data_c13z_r2_15n);
		PREMultiSet pre_data_c13z(singleset_vec_c13z, cam_, 13, "Nitroxide");
		pre_data_c13z.set_protein_mass(16.1);
		pre_data_c13z.set_temperature(298.0);
		pre_data_c13z.set_tau_c_limits(12.0e-9, 18.0e-9);

		pre_data_c13z.update_spin_coordinates(cam_);
		core::Real pre_score = pre_data_c13z.compute_pre_score_from_point_vector(radical_atom_wts_coords_);
		TR.Debug << "PRE score for PREMultiSet at spinlabel site " << pre_data_c13z.get_spinlabel_site_rsd() << " is " << pre_score << std::endl;
		pre_data_c13z.show(TR.Debug);
		TS_ASSERT_DELTA(pre_score,1.85747e-07,1.0e-3);
	}

	/// @brief test PRE score calculation with grid search
	void test_pre_score_calculation_with_gridsearch() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		TR << "Testing PREMultiSet score calculation with grid search" << std::endl;

		// Center grid search in MN ion
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("MN", 44), mid1_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("MN", 44), mid1_) );
		NMRGridSearchOP gs1( new NMRGridSearch(grid_atom1, grid_atom2, mid1_, 0, 0.5, 0, 1) );
		utility::vector1<PRESingleSetOP> singleset_vec_mn;
		PRESingleSetOP pres_mid1_mn_r1_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r1.txt", mid1_, 1.0, "R1", "CONST"));
		PRESingleSetOP pres_mid1_mn_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r2.txt", mid1_, 1.0, "R2", "CONST"));
		singleset_vec_mn.push_back(pres_mid1_mn_r1_1h);
		singleset_vec_mn.push_back(pres_mid1_mn_r2_1h);
		PREMultiSet pre_data_mid1_mn(singleset_vec_mn, mid1_, 44, "Mn2+", gs1);
		pre_data_mid1_mn.set_temperature(298.0);
		pre_data_mid1_mn.set_protein_mass(10.2);

		// Stepsize > out radius of gridsearch -> no traversal of metal coordinates
		NMRGridSearchOP gs2( new NMRGridSearch(grid_atom1, grid_atom2, mid1_, 0, 1.0, 0, 0.1));
		utility::vector1<PRESingleSetOP> singleset_vec_cu;
		PRESingleSetOP pres_mid1_cu_r1_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_cu2_1h_r1.txt", mid1_, 1.0, "R1", "CONST"));
		PRESingleSetOP pres_mid1_cu_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_cu2_1h_r2.txt", mid1_, 1.0, "R2", "CONST"));
		singleset_vec_cu.push_back(pres_mid1_cu_r1_1h);
		singleset_vec_cu.push_back(pres_mid1_cu_r2_1h);
		PREMultiSet pre_data_mid1_cu(singleset_vec_cu, mid1_, 44, "Cu2+", gs2);
		pre_data_mid1_cu.set_temperature(298.0);
		pre_data_mid1_cu.set_protein_mass(10.2);

		// First calculate the score for a fixed point which is already the metal ion position
		pre_data_mid1_mn.update_spin_coordinates(mid1_);
		gs1->set_grid_search_center(mid1_);
		core::Vector mn(gs1->get_grid_search_center());
		core::Real score_mn = pre_data_mid1_mn.compute_pre_score_from_single_point(mn);
		TR.Debug << "PRE score for Mn2+ ion at site " << pre_data_mid1_mn.get_spinlabel_site_rsd() << " is " << score_mn << std::endl;
		pre_data_mid1_mn.show(TR.Debug);

		pre_data_mid1_cu.update_spin_coordinates(mid1_);
		gs2->set_grid_search_center(mid1_);
		core::Vector cu(gs2->get_grid_search_center());
		core::Real score_cu = pre_data_mid1_cu.compute_pre_score_from_single_point(cu);
		TR.Debug << "PRE score for Cu2+ ion at site " << pre_data_mid1_cu.get_spinlabel_site_rsd() << " is " << score_cu << std::endl;
		pre_data_mid1_cu.show(TR.Debug);

		// Now, we are slightly off the center and optimize the position
		mn += core::Vector(0.5, 0.5, -0.5);
		core::Real score_mn_opt = pre_data_mid1_mn.compute_pre_score_from_single_point(mn, true);
		TS_ASSERT_DELTA(mn.x(), gs1->get_grid_search_center().x(), 1.0e-1);
		TS_ASSERT_DELTA(mn.y(), gs1->get_grid_search_center().y(), 1.0e-1);
		TS_ASSERT_DELTA(mn.z(), gs1->get_grid_search_center().z(), 1.0e-1);
		TR.Debug << "PRE score for Mn2+ ion at site " << pre_data_mid1_mn.get_spinlabel_site_rsd() << " is " << score_mn_opt << std::endl;
		pre_data_mid1_mn.show(TR.Debug);

		cu += core::Vector(0.5, -0.5, 0.5);
		core::Real score_cu_opt = pre_data_mid1_cu.compute_pre_score_from_single_point(cu, true);
		TS_ASSERT_DELTA(cu.x(), gs2->get_grid_search_center().x(), 1.0e-1);
		TS_ASSERT_DELTA(cu.y(), gs2->get_grid_search_center().y(), 1.0e-1);
		TS_ASSERT_DELTA(cu.z(), gs2->get_grid_search_center().z(), 1.0e-1);
		TR.Debug << "PRE score for Cu2+ ion at site " << pre_data_mid1_cu.get_spinlabel_site_rsd() << " is " << score_cu_opt << std::endl;
		pre_data_mid1_cu.show(TR.Debug);
	}

	/// @brief calculate the PRE score by building a spinlabel ensemble and by finding the metal
	///        position that best fits the PRE data by using a grid search
	void test_find_spinlabel_position_and_pre_score_calculation() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true -nmr:pre:nls_repeats 5");

		TR << "Testing optimization of spinlabel position and PREMultiSet score calculation" << std::endl;

		NMRSpinlabelOP spinlabel( new NMRSpinlabel("fa_standard", "R1A") );
		PRESingleSetOP pre_data_c13z_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_1h.txt", cam_, 1.0, "R2", "CONST"));
		PRESingleSetOP pre_data_c13z_r2_15n( new PRESingleSet("core/scoring/nmr/pre/sim_pre_c13z_r2_15n.txt", cam_, 1.0, "R2", "CONST"));
		utility::vector1<PRESingleSetOP> singleset_vec_c13z;
		singleset_vec_c13z.push_back(pre_data_c13z_r2_1h);
		singleset_vec_c13z.push_back(pre_data_c13z_r2_15n);
		PREMultiSet pre_data_c13z(singleset_vec_c13z, cam_, 13, "Nitroxide", spinlabel);
		pre_data_c13z.set_temperature(298.0);
		pre_data_c13z.set_protein_mass(16.1);

		// Center grid search in MN ion and make a grid of +/-2 Ang. around the metal ion
		core::id::AtomID grid_atom1( named_atom_id_to_atom_id(core::id::NamedAtomID("MN", 44), mid1_) );
		core::id::AtomID grid_atom2( named_atom_id_to_atom_id(core::id::NamedAtomID("MN", 44), mid1_) );
		NMRGridSearchOP gs1( new NMRGridSearch(grid_atom1, grid_atom2, mid1_, 0, 1.0, 0, 2.0) );
		utility::vector1<PRESingleSetOP> singleset_vec_mn;
		PRESingleSetOP pres_mid1_mn_r1_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r1.txt", mid1_, 1.0, "R1", "CONST"));
		PRESingleSetOP pres_mid1_mn_r2_1h( new PRESingleSet("core/scoring/nmr/pre/sim_pres_mn2_1h_r2.txt", mid1_, 1.0, "R2", "CONST"));
		singleset_vec_mn.push_back(pres_mid1_mn_r1_1h);
		singleset_vec_mn.push_back(pres_mid1_mn_r2_1h);
		PREMultiSet pre_data_mid1_mn(singleset_vec_mn, mid1_, 44, "Mn2+", gs1);
		pre_data_mid1_mn.set_temperature(298.0);
		pre_data_mid1_mn.set_protein_mass(10.2);

		WeightCoordVector r1a_pos;
		WeightCoordVector mn_pos;

		core::Real score_r1a = pre_data_c13z.find_para_ion_position_and_compute_pre_score(cam_, r1a_pos);
		TR.Debug << "PRE score for spinlabel " << pre_data_c13z.get_spinlabel()->get_code() << " at site " << pre_data_c13z.get_spinlabel_site_rsd()
			<< " is " << score_r1a << std::endl;
		TS_ASSERT_EQUALS(r1a_pos.size(), 19);
		pre_data_c13z.show(TR.Debug);

		core::Real score_mn = pre_data_mid1_mn.find_para_ion_position_and_compute_pre_score(mid1_, mn_pos);
		TR.Debug << "PRE score for Mn2+ ion at site " << pre_data_mid1_mn.get_spinlabel_site_rsd() << " is " << score_mn << std::endl;
		TS_ASSERT_EQUALS(mn_pos.size(), 1);
		TS_ASSERT_DELTA(mn_pos[1].second.x(), -11.750, 1.0e-1);
		TS_ASSERT_DELTA(mn_pos[1].second.y(),  11.145, 1.0e-1);
		TS_ASSERT_DELTA(mn_pos[1].second.z(),  -4.037, 1.0e-1);
		pre_data_mid1_mn.show(TR.Debug);
	}

};
