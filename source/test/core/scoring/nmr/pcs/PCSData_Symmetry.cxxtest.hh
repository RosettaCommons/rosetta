// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSData_Symmetry.cxxtest.hh
/// @brief   unit test for class PCSData that stores and handles all PCS data for all tagging sites and all lanthanides
/// @details Last modified: 10/03/16
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
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

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

static basic::Tracer TR("core.scoring.nmr.pcs.PCSData_Symmetry.cxxtest");

class PCSDataSymmetryTests : public CxxTest::TestSuite {

private:
	core::pose::Pose il10_;

public:
	/// @brief Setup Test
	void setUp() {
		using namespace core::conformation::symmetry;

		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(il10_, "core/scoring/nmr/pcs/2ilk_chainA.pdb", core::import_pose::PDB_file);
		// Make a symmetric pose
		SymmDataOP symmdata = SymmDataOP( new SymmData );
		symmdata->read_symmetry_data_from_file( "core/scoring/nmr/pcs/2ilk_chainAB.symm" );
		core::pose::symmetry::make_symmetric_pose( il10_, *symmdata );
		runtime_assert( core::pose::symmetry::is_symmetric( il10_ ) );
	}

	void tearDown() {
		il10_.clear();
	}

	/// @brief PCS score calculation for symmetric protein by calling PCSData's compute_score_all_tags() method.
	///        The symmetric spins get deduced automatically from the symmetry info object
	///        and only the spins of the ASU are provided in the input file.
	void test_PCSData_score_calculation_all_tags_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr::pcs;

		// Important: use this option to turn automatic symmetry handling on
		// Sets also fixed RG seed
		core_init_with_additional_options("-nmr:pcs:use_symmetry_calc");

		PCSData pcs_data("core/scoring/nmr/pcs/pcs_data_2ilk_asu_input_file.txt", il10_);
		utility::vector1< utility::vector1< PCSTensorCOP > > tensors_all_tags_and_lanthanides(pcs_data.get_number_tags());
		utility::vector1< core::Real > scores_all_tags(pcs_data.get_number_tags());
		core::Real total_score = pcs_data.compute_score_all_tags( il10_, scores_all_tags, tensors_all_tags_and_lanthanides);

		// Show and check the results (For brevity, only the metal coordinates)
		TR.Debug << "Total PCS score is " << total_score << std::endl;
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().z(), 58.078, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().z(), 58.078, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().z(), 58.078, 1e-1);
	}

	/// @brief PCS score calculation for symmetric protein by calling PCSData's compute_score_all_tags() method.
	///        Here, the symmetric spins are explicitly declared in the input file
	///        and the PCS is calculated by sum averaging.
	void test_PCSData_score_calculation_all_tags_manual_symmetry_handling() {
		using namespace core::scoring::nmr::pcs;

		// Set fixed RG seed for NLS fitting
		initialize_rng();

		PCSData pcs_data("core/scoring/nmr/pcs/pcs_data_2ilk_2su_input_file.txt", il10_);
		// In the input file, sum averaging of the PCS values of symmetric spins is turned on.
		// We test here that this is indeed true.
		utility::vector1< PCSMultiSetOP > & pcs_multiset_vec = pcs_data.get_pcs_multiset_vec();
		core::Size number_tags = pcs_data.get_number_tags();
		for ( core::Size i = 1; i <= number_tags; ++i ) {
			utility::vector1< PCSSingleSetOP > & pcs_singleset_vec = pcs_multiset_vec[i]->get_pcs_singleset_vec();
			core::Size number_lanthanides = pcs_multiset_vec[i]->get_number_metal_ions();
			for ( core::Size j = 1; j <= number_lanthanides; ++j ) {
				TS_ASSERT(pcs_singleset_vec[j]->get_averaging_type() == core::scoring::nmr::SUM);
			}
		}

		utility::vector1< utility::vector1< PCSTensorCOP > > tensors_all_tags_and_lanthanides(number_tags);
		utility::vector1< core::Real > scores_all_tags(number_tags);
		core::Real total_score = pcs_data.compute_score_all_tags( il10_, scores_all_tags, tensors_all_tags_and_lanthanides);

		// Show and check the results (For brevity, only the metal coordinates)
		TR.Debug << "Total PCS score is " << total_score << std::endl;
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][1]->get_metal_center().z(), 58.078, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][2]->get_metal_center().z(), 58.078, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().x(),  2.776, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().y(), 41.602, 1e-1);
		TS_ASSERT_DELTA(tensors_all_tags_and_lanthanides[1][3]->get_metal_center().z(), 58.078, 1e-1);

	}
};
