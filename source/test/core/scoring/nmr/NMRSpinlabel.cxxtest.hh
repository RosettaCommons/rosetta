// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRSpinlabel.cxxtest.hh
/// @brief   unit test for class NMRSpinlabel.
/// @details Last modified: 10/05/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>

static basic::Tracer TR("core.scoring.nmr.NMRSpinlabel.cxxtest");

class NMRSpinlabelTests : public CxxTest::TestSuite {

private:
	core::pose::Pose cam_;
	std::string ensemble_conformers_database_file_;

public:
	/// @brief Setup Test
	void setUp() {

		// Initialize core & options system
		core_init();

		// Load pose from pdb
		core::import_pose::pose_from_file(cam_, "core/scoring/nmr/1cdl.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		cam_.clear();
	}

	/// @brief test creation of NMRSpinlabel object
	///        and access of NMRSpinlabel data
	void test_creation_and_access() {
		using namespace core::scoring::nmr;

		TR << "Testing instantiation of NMRSpinlabel." << std::endl;

		// Test creation of NMRSpinlabel object and some data members that are automatically set
		NMRSpinlabelOP spl( new NMRSpinlabel("fa_standard", "R1A") );
		TS_ASSERT_EQUALS(spl->get_name(), "R1A");
		TS_ASSERT_EQUALS(spl->get_code(),"R1A");
		TS_ASSERT_EQUALS(spl->get_radical_atom(), "O1");
		TS_ASSERT_EQUALS(spl->get_max_ensemble_size(), 20);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), 54);
		TS_ASSERT_EQUALS(spl->get_highres_conformer_filter_type(), NMRSpinlabel::BUMP_ENERGY);
		TS_ASSERT_DELTA(spl->get_boltzman_factor(), 2.0, 1.0e-6);

		// Test that NMRDummySpinlabelConformers have been aligned to the coordinate origin
		NMRDummySpinlabelEnsembleOP ndsl_ensemble(spl->get_dummy_ensemble());
		NMRDummySpinlabelConformerOP ndsl_conformer1(ndsl_ensemble->get_conformer_table()[1]);
		NMRDummySpinlabelConformerOP ndsl_conformer2(ndsl_ensemble->get_conformer_table()[2]);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CA").get_coordinates().x(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CA").get_coordinates().y(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CA").get_coordinates().z(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CB").get_coordinates().x(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CB").get_coordinates().y(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("CB").get_coordinates().z(),1.52871, 1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("C").get_coordinates().x(), 1.42489, 1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("C").get_coordinates().y(), 0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer1->get_atom_table().at("C").get_coordinates().z(),-0.53602, 1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CA").get_coordinates().x(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CA").get_coordinates().y(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CA").get_coordinates().z(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CB").get_coordinates().x(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CB").get_coordinates().y(),0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("CB").get_coordinates().z(),1.52871, 1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("C").get_coordinates().x(), 1.42489, 1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("C").get_coordinates().y(), 0.0,     1.0e-2);
		TS_ASSERT_DELTA(ndsl_conformer2->get_atom_table().at("C").get_coordinates().z(),-0.53602, 1.0e-2);

		// Test that RMSD matrix between NMRDummySpinlabelConformers has been created correctly
		utility::vector1<utility::vector1<core::Real>> const & rmsd_mat = ndsl_ensemble->get_rmsd_mat();
		TS_ASSERT_DELTA(rmsd_mat[1][2], 1.93813, 1.0e-2);
		TS_ASSERT_DELTA(rmsd_mat[1][3], 3.80159, 1.0e-2);
		TS_ASSERT_DELTA(rmsd_mat[1][4], 3.67895, 1.0e-2);
		core::conformation::ResidueCOP rsd1 = ndsl_conformer1->get_residue();
		core::conformation::ResidueCOP rsd2 = ndsl_conformer2->get_residue();
		core::Real rmsd12 = core::scoring::automorphic_rmsd(*rsd1,*rsd2,false);
		TS_ASSERT_DELTA(rmsd12, 1.93813, 1.0e-2);
	}

	/// @brief test filtering NMRDummySpinlabelEnsemble
	///        by neighbor count
	void test_filter_spinlabel_ensemble_by_distance_check() {
		using namespace core::scoring::nmr;
		using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

		core_init_with_additional_options("-nmr:spinlabel:max_ensemble_size 100");
		TR << "Testing NMRSpinlabel neighbor count filter." << std::endl;

		NMRSpinlabelOP spl( new NMRSpinlabel("fa_standard", "R1A") );
		core::Size n_good_conformers(19);
		core::Size n_conformers(54);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_conformers);
		WeightCoordVector radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_distance_check(cam_, 13);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);

		NMRDummySpinlabelEnsembleOP ndsl_ensemble(spl->get_dummy_ensemble());
		utility::vector1<NMRDummySpinlabelConformerOP> & ndsl_conformers = ndsl_ensemble->get_conformer_table();
		std::list<core::Size> ndsl_good_conformers({3, 6, 13, 14, 17, 18, 19, 20, 21, 24, 26, 30, 32, 34, 36, 38, 40, 45, 50});
		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( std::find(ndsl_good_conformers.begin(), ndsl_good_conformers.end(), conformer_id) != ndsl_good_conformers.end() ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

		// Also test ensemble creation at N- and C-terminus
		radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_distance_check(cam_, 1);
		n_good_conformers=32;
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);
		ndsl_good_conformers.clear();
		ndsl_good_conformers = std::list<core::Size>({7, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19, 20, 21, 23, 25, 27, 29, 30, 31, 32, 33, 36, 38, 39, 41, 42, 43, 44, 45, 47, 49, 50});
		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( std::find(ndsl_good_conformers.begin(), ndsl_good_conformers.end(), conformer_id) != ndsl_good_conformers.end() ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

		radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_distance_check(cam_, cam_.total_residue());
		n_good_conformers=31;
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);
		ndsl_good_conformers.clear();
		ndsl_good_conformers = std::list<core::Size>({1, 2, 3, 6, 7, 9, 11, 13, 14, 15, 18, 20, 22, 24, 26, 27, 28, 30, 32, 34, 40, 41, 42, 43, 44, 46, 48, 51, 52, 53, 54});
		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( std::find(ndsl_good_conformers.begin(), ndsl_good_conformers.end(), conformer_id) != ndsl_good_conformers.end() ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

	}

	/// @brief test filtering NMRDummySpinlabelEnsemble
	///        by neighbor count and additional clustering
	void test_filter_spinlabel_ensemble_by_distance_check_and_clustering() {
		using namespace core::scoring::nmr;
		using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

		core_init_with_additional_options("-nmr:spinlabel:max_ensemble_size 10");
		TR << "Testing NMRSpinlabel neighbor count filter and clustering." << std::endl;

		NMRSpinlabelOP spl( new NMRSpinlabel("fa_standard", "R1A") );
		core::Size n_good_conformers(10);
		core::Size n_conformers(54);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_conformers);
		WeightCoordVector radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_distance_check(cam_, 13);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);

	}

	/// @brief test filtering NMRDummySpinlabelEnsemble
	///        by mask
	void test_filter_spinlabel_ensemble_by_mask() {
		using namespace core::scoring::nmr;
		using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

		core_init_with_additional_options("-nmr:spinlabel:max_ensemble_size 100");

		TR << "Testing NMRSpinlabel filter by masking" << std::endl;

		NMRSpinlabelOP spl( new NMRSpinlabel("fa_standard", "R1A") );
		core::Size n_good_conformers(19);
		core::Size n_conformers(54);
		std::list<core::Size> ndsl_good_conformers({3, 8, 10, 13, 14, 17, 19, 20, 21, 27, 30, 32, 34, 36, 38, 40, 45, 46, 50});
		utility::vector1<bool> mask(n_conformers,false);
		for ( auto i : ndsl_good_conformers ) { mask[i]=true; }

		utility::vector1<core::Real> conformer_energies(n_conformers);
		conformer_energies[1]=1709.6;     conformer_energies[2]=1848.93;   conformer_energies[3]=-4.48112;
		conformer_energies[4]=2676.86;    conformer_energies[5]=2753.81;   conformer_energies[6]=62.5831;
		conformer_energies[7]=4145.84;    conformer_energies[8]=-4.09032;  conformer_energies[9]=4440.55;
		conformer_energies[10]=-0.110651; conformer_energies[11]=1098.08;  conformer_energies[12]=1297.06;
		conformer_energies[13]=-1.74973;  conformer_energies[14]=-1.7598;  conformer_energies[15]=201.854;
		conformer_energies[16]=1173.04;   conformer_energies[17]=-2.24815; conformer_energies[18]=186.59;
		conformer_energies[19]=-3.82332;  conformer_energies[20]=2.17002;  conformer_energies[21]=-2.32267;
		conformer_energies[22]=2058.16;   conformer_energies[23]=1628.26;  conformer_energies[24]=26.6825;
		conformer_energies[25]=1473.36;   conformer_energies[26]=8.58285;  conformer_energies[27]=-2.35253;
		conformer_energies[28]=1555.81;   conformer_energies[29]=1904.24;  conformer_energies[30]=-4.85482;
		conformer_energies[31]=2755.3;    conformer_energies[32]=-4.79701; conformer_energies[33]=4559.98;
		conformer_energies[34]=-5.16737;  conformer_energies[35]=5050.94;  conformer_energies[36]=-5.55905;
		conformer_energies[37]=5239.14;   conformer_energies[38]=-5.61152; conformer_energies[39]=3076.66;
		conformer_energies[40]=-5.32601;  conformer_energies[41]=902.986;  conformer_energies[42]=4205.99;
		conformer_energies[43]=472.967;   conformer_energies[44]=3631.05;  conformer_energies[45]=-5.15183;
		conformer_energies[46]=-6.72808;  conformer_energies[47]=2927.91;  conformer_energies[48]=10.9597;
		conformer_energies[49]=1784.22;   conformer_energies[50]=-5.4108;  conformer_energies[51]=42.3042;
		conformer_energies[52]=21.1767;   conformer_energies[53]=17.189;   conformer_energies[54]=468.664;

		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_conformers);
		WeightCoordVector radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_mask(cam_, 13, mask, conformer_energies);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);

		NMRDummySpinlabelEnsembleOP ndsl_ensemble(spl->get_dummy_ensemble());
		utility::vector1<NMRDummySpinlabelConformerOP> & ndsl_conformers = ndsl_ensemble->get_conformer_table();
		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( mask[ conformer_id ] ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

		// Also test ensemble creation at N- and C-terminus
		n_good_conformers = 30;
		std::fill(mask.begin(), mask.end(), false);
		ndsl_good_conformers.clear();
		ndsl_good_conformers = std::list<core::Size>({8,10,12,13,14,16,17,18,19,20,21,23,25,26,27,29,30,31,32,33,36,38,39,41,42,44,45,47,49,50});
		for ( auto i : ndsl_good_conformers ) { mask[i]=true; }
		conformer_energies[1]=2912.82;     conformer_energies[2]=3251.47;     conformer_energies[3]=171.863;
		conformer_energies[4]=447.198;     conformer_energies[5]=1249.27;     conformer_energies[6]=223.699;
		conformer_energies[7]=632.666;     conformer_energies[8]=-0.360856;   conformer_energies[9]=424.495;
		conformer_energies[10]=-0.466888;  conformer_energies[11]=2776.26;    conformer_energies[12]=-0.413284;
		conformer_energies[13]=-0.531125;  conformer_energies[14]=-0.730523;  conformer_energies[15]=2975.61;
		conformer_energies[16]=-0.420579;  conformer_energies[17]=-1.87701;   conformer_energies[18]=-3.24149;
		conformer_energies[19]=-1.87912;   conformer_energies[20]=-3.01552;   conformer_energies[21]=-2.28207;
		conformer_energies[22]=898.854;    conformer_energies[23]=-1.89915;   conformer_energies[24]=458.665;
		conformer_energies[25]=-1.88966;   conformer_energies[26]=-0.0749706; conformer_energies[27]=-2.22106;
		conformer_energies[28]=635.705;    conformer_energies[29]=-0.759497;  conformer_energies[30]=-0.766275;
		conformer_energies[31]=-1.09892;   conformer_energies[32]=-0.728941;  conformer_energies[33]=-1.45215;
		conformer_energies[34]=64.8607;    conformer_energies[35]=207.464;    conformer_energies[36]=-0.834449;
		conformer_energies[37]=811.657;    conformer_energies[38]=-0.74937;   conformer_energies[39]=-0.844671;
		conformer_energies[40]=104.054;    conformer_energies[41]=-1.90273;   conformer_energies[42]=-1.7866;
		conformer_energies[43]=86.7842;    conformer_energies[44]=-1.22912;   conformer_energies[45]=-1.42966;
		conformer_energies[46]=1302.75;    conformer_energies[47]=-1.12692;   conformer_energies[48]=1126.36;
		conformer_energies[49]=-1.12267;   conformer_energies[50]=-1.39062;   conformer_energies[51]=1772.98;
		conformer_energies[52]=2243.21;    conformer_energies[53]=2529.8;     conformer_energies[54]=1845.13;

		radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_mask(cam_, 1, mask, conformer_energies);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);

		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( mask[ conformer_id ] ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

		n_good_conformers = 28;
		std::fill(mask.begin(), mask.end(), false);
		ndsl_good_conformers.clear();
		ndsl_good_conformers = std::list<core::Size>({3,6,7,9,11,13,14,15,17,18,19,20,21,24,26,27,30,32,34,40,41,43,46,48,51,52,53,54});
		for ( auto i : ndsl_good_conformers ) { mask[i]=true; }
		conformer_energies[1]=139.392;    conformer_energies[2]=587.837;    conformer_energies[3]=-2.09657;
		conformer_energies[4]=2406.33;    conformer_energies[5]=2168.21;    conformer_energies[6]=-1.62763;
		conformer_energies[7]=-2.74392;   conformer_energies[8]=160.751;    conformer_energies[9]=-2.93765;
		conformer_energies[10]=600.191;   conformer_energies[11]=3.71585;   conformer_energies[12]=525.675;
		conformer_energies[13]=-2.17413;  conformer_energies[14]=-1.47354;  conformer_energies[15]=-1.53328;
		conformer_energies[16]=398.883;   conformer_energies[17]=-3.64353;  conformer_energies[18]=-0.876046;
		conformer_energies[19]=-2.5676;   conformer_energies[20]=-1.94987;  conformer_energies[21]=-1.84675;
		conformer_energies[22]=1636.71;   conformer_energies[23]=411.828;   conformer_energies[24]=-1.78848;
		conformer_energies[25]=568.863;   conformer_energies[26]=-1.33285;  conformer_energies[27]=-2.60089;
		conformer_energies[28]=1521.67;   conformer_energies[29]=1138.79;   conformer_energies[30]=-2.94093;
		conformer_energies[31]=799.798;   conformer_energies[32]=-3.59993;  conformer_energies[33]=4234.98;
		conformer_energies[34]=-2.66184;  conformer_energies[35]=4929.14;   conformer_energies[36]=360.53;
		conformer_energies[37]=4019;      conformer_energies[38]=187.581;   conformer_energies[39]=2320.43;
		conformer_energies[40]=-2.76716;  conformer_energies[41]=-2.23385;  conformer_energies[42]=23.1074;
		conformer_energies[43]=-1.99856;  conformer_energies[44]=23.5578;   conformer_energies[45]=336.487;
		conformer_energies[46]=-2.9466;   conformer_energies[47]=2317.9;    conformer_energies[48]=-2.97839;
		conformer_energies[49]=1956.79;   conformer_energies[50]=112.562;   conformer_energies[51]=-1.12033;
		conformer_energies[52]=-0.96885;  conformer_energies[53]=-0.940899; conformer_energies[54]=-1.60185;

		radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_mask(cam_, 1, mask, conformer_energies);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_good_conformers);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_good_conformers);

		for ( core::Size conformer_id(1); conformer_id<=n_conformers; conformer_id++ ) {
			if ( mask[ conformer_id ] ) {
				TS_ASSERT( !ndsl_conformers[ conformer_id ]->has_clash() );
			} else {
				TS_ASSERT( ndsl_conformers[ conformer_id ]->has_clash() );
			}
		}

	}

	/// @brief test filtering NMRDummySpinlabelEnsemble
	///        by mask and clustering
	void test_filter_spinlabel_ensemble_by_mask_and_clustering() {
		using namespace core::scoring::nmr;
		using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

		core_init_with_additional_options("-nmr:spinlabel:max_ensemble_size 10");

		TR << "Testing NMRSpinlabel filter by masking and clustering" << std::endl;

		NMRSpinlabelOP spl( new NMRSpinlabel("fa_standard", "R1A") );
		core::Size n_good_conformers(19);
		core::Size n_clusters(10);
		core::Size n_conformers(54);
		int conformer_ids[] = {3, 8, 10, 13, 14, 17, 19, 20, 21, 27, 30, 32, 34, 36, 38, 40, 45, 46, 50};
		utility::vector1<bool> mask(n_conformers,false);
		for ( core::Size i(0); i<n_good_conformers; i++ ) { mask[ conformer_ids[i] ]=true; }

		utility::vector1<core::Real> conformer_energies(n_conformers);
		conformer_energies[1]=1709.6;     conformer_energies[2]=1848.93;   conformer_energies[3]=-4.48112;
		conformer_energies[4]=2676.86;    conformer_energies[5]=2753.81;   conformer_energies[6]=62.5831;
		conformer_energies[7]=4145.84;    conformer_energies[8]=-4.09032;  conformer_energies[9]=4440.55;
		conformer_energies[10]=-0.110651; conformer_energies[11]=1098.08;  conformer_energies[12]=1297.06;
		conformer_energies[13]=-1.74973;  conformer_energies[14]=-1.7598;  conformer_energies[15]=201.854;
		conformer_energies[16]=1173.04;   conformer_energies[17]=-2.24815; conformer_energies[18]=186.59;
		conformer_energies[19]=-3.82332;  conformer_energies[20]=2.17002;  conformer_energies[21]=-2.32267;
		conformer_energies[22]=2058.16;   conformer_energies[23]=1628.26;  conformer_energies[24]=26.6825;
		conformer_energies[25]=1473.36;   conformer_energies[26]=8.58285;  conformer_energies[27]=-2.35253;
		conformer_energies[28]=1555.81;   conformer_energies[29]=1904.24;  conformer_energies[30]=-4.85482;
		conformer_energies[31]=2755.3;    conformer_energies[32]=-4.79701; conformer_energies[33]=4559.98;
		conformer_energies[34]=-5.16737;  conformer_energies[35]=5050.94;  conformer_energies[36]=-5.55905;
		conformer_energies[37]=5239.14;   conformer_energies[38]=-5.61152; conformer_energies[39]=3076.66;
		conformer_energies[40]=-5.32601;  conformer_energies[41]=902.986;  conformer_energies[42]=4205.99;
		conformer_energies[43]=472.967;   conformer_energies[44]=3631.05;  conformer_energies[45]=-5.15183;
		conformer_energies[46]=-6.72808;  conformer_energies[47]=2927.91;  conformer_energies[48]=10.9597;
		conformer_energies[49]=1784.22;   conformer_energies[50]=-5.4108;  conformer_energies[51]=42.3042;
		conformer_energies[52]=21.1767;   conformer_energies[53]=17.189;   conformer_energies[54]=468.664;
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_conformers);
		WeightCoordVector radical_atom_wghts_coords = spl->filter_spinlabel_ensemble_by_mask(cam_, 13, mask, conformer_energies);
		TS_ASSERT_EQUALS(spl->get_current_ensemble_size(), n_clusters);
		TS_ASSERT_EQUALS(radical_atom_wghts_coords.size(), n_clusters);
	}
};
