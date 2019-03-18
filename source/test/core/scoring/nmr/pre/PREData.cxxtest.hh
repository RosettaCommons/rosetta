// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREData.cxxtest.hh
/// @brief   unit test for class PREData that stores and handles data of multiple PRE experiments collected for multiple different spinlabels
///          those could be either different spinlabel tagging positions or chemically different spinlabels
/// @details Last modified: 10/16/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/pre/PREData.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/io/nmr/ParaIon.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

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

static basic::Tracer TR("core.scoring.nmr.pre.PREData.cxxtest");

class PREDataTests : public CxxTest::TestSuite {

	using WeightCoordVector = core::scoring::nmr::pre::PREMultiSet::WeightCoordVector;

private:
	core::pose::Pose cam_;
	core::pose::Pose mid1_;
	utility::vector1< WeightCoordVector > radical_atm_wts_coos_all_sites_;

public:

	/// @brief Setup Test
	void setUp() {

		// Initialize core & options system
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		// Load pose from pdb
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
		core::import_pose::pose_from_file(mid1_, "core/scoring/nmr/pre/3v1e_AB_Mn.pdb", core::import_pose::PDB_file);

		radical_atm_wts_coos_all_sites_.resize(4);
		// C13Z
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(81.299, 33.311,  7.942)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(92.649, 31.000,  8.182)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(91.389, 32.009, 10.363)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.609, 34.451,  9.000)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.679, 33.309, 10.941)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(87.027, 35.888,  7.247)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(87.214, 35.445,  4.572)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.878, 34.328,  7.101)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.093, 32.216, 10.080)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.145, 32.886, 10.977)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(83.631, 34.837,  7.928)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.025, 35.951,  7.301)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(84.344, 30.100, 10.762)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.231, 32.869,  9.176)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.496, 33.889,  7.204)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.633, 30.224,  9.316)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.796, 34.025,  8.689)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(84.720, 29.332, 11.161)));
		radical_atm_wts_coos_all_sites_[1].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(90.801, 33.438,  7.847)));
		// D46Z
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(84.646,  2.609, -4.577)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(85.103,  4.110, -6.193)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(95.846, -0.192, -3.610)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(94.602, -0.249, -6.021)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(92.103, -0.465,  0.515)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(89.100, -0.552, -5.461)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(89.908,  0.702, -7.149)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.636, -2.255, -3.978)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(88.612, -2.307, -1.260)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(83.729,  1.768, -3.969)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(93.311,  0.086, -5.869)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(91.293,  0.317,  1.930)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(85.852,  5.013, -7.002)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(91.407,  0.452, -7.006)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(86.112,  0.240, -4.685)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(87.716, -1.871, -4.110)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(93.038, -0.679, -5.094)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(92.578, -1.974, -3.329)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(95.799,  1.492, -6.577)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(95.237,  4.037, -7.327)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(91.201, -1.158, -4.927)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(97.011,  0.851,  1.667)));
		radical_atm_wts_coos_all_sites_[2].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(93.128, -1.615, -3.857)));
		// E83Z
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(81.122, 17.926, -3.046)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.735, 12.985, -5.051)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.115, 13.104, -2.890)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(74.639,  4.989, -7.233)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(76.413,  4.655, -5.206)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.897,  7.833, -5.886)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.132,  7.547, -3.787)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(81.662,  7.519, -8.039)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.888,  8.633,-10.395)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(83.714, 13.173, -6.014)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(77.261,  5.696, -5.173)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(80.701, 15.138, -1.788)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(81.304, 12.973, -1.725)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(78.984,  6.571, -4.026)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.724, 10.430, -6.168)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(82.268,  8.263, -7.735)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(78.727, 10.966, -2.333)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(77.624,  5.637, -6.234)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(78.101,  5.727, -8.417)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(79.601, 12.832, -3.222)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(74.750,  4.962, -3.822)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(74.193,  6.764, -1.873)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(79.312,  6.424, -6.639)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(77.929, 10.871, -1.837)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(75.790, 11.399, -1.375)));
		radical_atm_wts_coos_all_sites_[3].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(77.634,  5.443, -7.776)));
		// E119Z
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(72.709, 37.240,  4.032)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(72.689, 37.275,  1.781)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(61.669, 40.727,  4.475)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(63.508, 42.113,  3.041)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(68.666, 40.959,  4.696)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(68.328, 41.287,  2.493)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(68.730, 41.158,  6.992)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(68.033, 39.419,  8.961)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(64.711, 41.530,  3.168)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(66.847, 41.653,  2.477)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(71.337, 39.321,  5.281)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(69.650, 40.793,  6.814)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(68.916, 37.919,  0.471)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(64.773, 41.537,  4.290)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(64.756, 41.253,  6.508)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(70.304, 36.441,  1.430)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(62.492, 41.416,  1.237)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(63.219, 39.931, -0.913)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(66.503, 41.443,  5.084)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(66.751, 36.516, -1.283)));
		radical_atm_wts_coos_all_sites_[4].push_back(std::make_pair(1, numeric::xyzVector< core::Real >(64.363, 41.435,  5.785)));
	}

	void tearDown() {
		cam_.clear();
		mid1_.clear();
		radical_atm_wts_coos_all_sites_.clear();
	}

	/// @brief PREData instantiation from input file using the PREData constructor
	void test_PREData_instantiation() {
		using namespace core::scoring::nmr::pre;

		TR << "Testing PREData creation" << std::endl;

		PREData pre_data_cam("core/scoring/nmr/pre/pre_data_1cdl_input.txt", cam_);
		pre_data_cam.show(TR.Debug);
		TS_ASSERT_EQUALS(pre_data_cam.get_number_spinlabel_sites(), 8);

		PREData pre_data_mid1("core/scoring/nmr/pre/pre_data_mid1_input.txt", mid1_);
		pre_data_mid1.show(TR.Debug);
		TS_ASSERT_EQUALS(pre_data_mid1.get_number_spinlabel_sites(), 2);

	}

	/// @brief PREData instantiation using the NMRFactory
	void test_PREData_creation_by_NMRFactory() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::pre;

		TR << "Testing PREData creation" << std::endl;

		PREDataOP pre_data_cam(utility::pointer::static_pointer_cast< PREData >( NMRDataFactory::get_instance()->get_nmr_data("PRE", "core/scoring/nmr/pre/pre_data_1cdl_input.txt", cam_) ) );
		pre_data_cam->show(TR.Debug);
		TS_ASSERT_EQUALS(pre_data_cam->get_number_spinlabel_sites(), 8);

		PREDataOP pre_data_mid1(utility::pointer::static_pointer_cast< PREData >( NMRDataFactory::get_instance()->get_nmr_data("PRE", "core/scoring/nmr/pre/pre_data_mid1_input.txt", mid1_) ) );
		pre_data_mid1->show(TR.Debug);
		TS_ASSERT_EQUALS(pre_data_mid1->get_number_spinlabel_sites(), 2);
	}

	/// @brief PRE score calculation with predefined spinlabel coordinates
	void test_pre_score_calculation_with_explicit_spinlabel_coordinates() {
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		TR << "Testing PREData score calculation with spinlabel" << std::endl;

		PREData pre_data("core/scoring/nmr/pre/pre_data_1cdl_input.txt", cam_);
		utility::vector1< PREMultiSetOP > & multiset_vec = pre_data.get_pre_multiset_vec();
		utility::vector1< core::Real > pre_scores_all_sl_sites(pre_data.get_number_spinlabel_sites());

		core::Size n_spinlabel_sites(radical_atm_wts_coos_all_sites_.size()); // 4 spin-label sites
		for ( core::Size i = 1; i <= n_spinlabel_sites; ++i ) {
			for ( core::Size j = 2*i-1, j_end = 2*i; j <= j_end; ++j ) { // every spin-label site has two PRE datasets
				multiset_vec[j]->update_spin_coordinates(cam_);
				pre_scores_all_sl_sites[j] = multiset_vec[j]->compute_pre_score_from_point_vector(radical_atm_wts_coos_all_sites_[i]);
				TR.Debug << "PRE score for spinlabel " << multiset_vec[j]->get_spinlabel()->get_code() << " at residue "
					<< multiset_vec[j]->get_spinlabel_site_rsd() << ": " << pre_scores_all_sl_sites[j] << std::endl;
				multiset_vec[j]->show(TR.Debug);
			}
		}

		//Since we don't use the spinlabel here but provide fixed input coordinates
		//the ensemble size of the spinlabel data member is still the original size
		TS_ASSERT_EQUALS(multiset_vec[1]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[2]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[3]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[4]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[5]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[6]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[7]->get_spinlabel()->get_current_ensemble_size(),54);
		TS_ASSERT_EQUALS(multiset_vec[8]->get_spinlabel()->get_current_ensemble_size(),54);
	}

	/// @brief PRE score calculation with implicit spinlabel method
	void test_pre_score_calculation_with_implicit_spinlabel_ensemble() {
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true -nmr:spinlabel:max_ensemble_size 60");

		TR << "Testing PREData score calculation with spinlabel" << std::endl;

		PREData pre_data("core/scoring/nmr/pre/pre_data_1cdl_input.txt", cam_);
		utility::vector1< PREMultiSetOP > & multiset_vec = pre_data.get_pre_multiset_vec();
		utility::vector1< core::Real > pre_scores_all_sl_sites(pre_data.get_number_spinlabel_sites());

		pre_data.compute_score_all_spinlabel(cam_, pre_scores_all_sl_sites);

		for ( core::Size i = 1; i <= pre_data.get_number_spinlabel_sites(); ++i ) {
			TR.Debug << "PRE score for spinlabel " << multiset_vec[i]->get_spinlabel()->get_code() << " at residue "
				<< multiset_vec[i]->get_spinlabel_site_rsd() << ": " << pre_scores_all_sl_sites[i] << std::endl;
			multiset_vec[i]->show(TR.Debug);
		}

		TS_ASSERT_EQUALS(multiset_vec[1]->get_spinlabel()->get_current_ensemble_size(),19);
		TS_ASSERT_EQUALS(multiset_vec[2]->get_spinlabel()->get_current_ensemble_size(),19);
		TS_ASSERT_EQUALS(multiset_vec[3]->get_spinlabel()->get_current_ensemble_size(),29);
		TS_ASSERT_EQUALS(multiset_vec[4]->get_spinlabel()->get_current_ensemble_size(),29);
		TS_ASSERT_EQUALS(multiset_vec[5]->get_spinlabel()->get_current_ensemble_size(),29);
		TS_ASSERT_EQUALS(multiset_vec[6]->get_spinlabel()->get_current_ensemble_size(),29);
		TS_ASSERT_EQUALS(multiset_vec[7]->get_spinlabel()->get_current_ensemble_size(),25);
		TS_ASSERT_EQUALS(multiset_vec[8]->get_spinlabel()->get_current_ensemble_size(),25);
	}

	/// @brief PRE score calculation for paramagnetic metal ion with gridsearch
	void test_pre_score_calculation_with_metal_ion_gridsearch() {
		using namespace core::scoring::nmr::pre;

		// Turn on normalization of input PREs by experiment StdDev; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true");

		TR << "Testing PREData score calculation with metal ion grid search" << std::endl;

		PREData pre_data("core/scoring/nmr/pre/pre_data_mid1_input.txt", mid1_);
		utility::vector1< PREMultiSetOP > & multiset_vec = pre_data.get_pre_multiset_vec();
		utility::vector1< core::Real > pre_scores_all_ions(pre_data.get_number_spinlabel_sites());

		pre_data.compute_score_all_spinlabel(mid1_, pre_scores_all_ions);

		for ( core::Size i = 1; i <= pre_data.get_number_spinlabel_sites(); ++i ) {
			TR.Debug << "PRE score for ion " << multiset_vec[i]->get_ion_type()->get_ion_label() << " at position "
				<< multiset_vec[i]->get_spinlabel_site_rsd() << ": " << pre_scores_all_ions[i] << std::endl;
			multiset_vec[i]->show(TR.Debug);
		}

		TS_ASSERT_DELTA(multiset_vec[1]->get_gridsearch_iterator()->get_best_grid_point().x(), -11.750, 1.0e-1);
		TS_ASSERT_DELTA(multiset_vec[1]->get_gridsearch_iterator()->get_best_grid_point().y(),  11.145, 1.0e-1);
		TS_ASSERT_DELTA(multiset_vec[1]->get_gridsearch_iterator()->get_best_grid_point().z(),  -4.037, 1.0e-1);
		TS_ASSERT_DELTA(multiset_vec[2]->get_gridsearch_iterator()->get_best_grid_point().x(), -11.750, 1.0e-1);
		TS_ASSERT_DELTA(multiset_vec[2]->get_gridsearch_iterator()->get_best_grid_point().y(),  11.145, 1.0e-1);
		TS_ASSERT_DELTA(multiset_vec[2]->get_gridsearch_iterator()->get_best_grid_point().z(),  -4.037, 1.0e-1);

	}

};
