// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCData.cxxtest.hh
/// @brief   unit test for class RDCData that stores and handles all rdc data for multiple alignment media
/// @details Last modified: 07/09/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/rdc/RDCData.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.hh>
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
#include <utility>

static basic::Tracer TR("core.scoring.nmr.rdc.RDCData_Symmetry.cxxtest");

class RDCDataSymmetryTests : public CxxTest::TestSuite {

private: // Data

	core::pose::Pose tolR_;

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

	}

	void tearDown() {
		tolR_.clear();
	}

	/// @brief Test RDC score calculation by calling RDCData's compute_score_all_media() method.
	///        The symmetric spins get deduced automatically from the symmetry info object
	///        and only the spins of the ASU are provided in the input file.
	void test_RDCData_score_calculation_all_media_automatic_symmetry_deduction() {
		using namespace core::scoring::nmr::rdc;

		// Important: use this option to turn automatic symmetry handling on
		// Sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH -nmr:rdc:use_symmetry_calc -nmr:rdc:multiset_weights 2.0 1.5");

		TR << "Testing RDCData score calculation for symmetric pose" << std::endl;

		// Setup RDC data
		RDCData rdc_data_all("core/scoring/nmr/rdc/rdc_data_2jwk_asu_input_file.txt", tolR_);

		utility::vector1< core::Real > scores_all_media;
		utility::vector1< RDCTensorCOP > tensors_all_media;
		core::Real total_score = rdc_data_all.compute_score_all_media(tolR_, scores_all_media, tensors_all_media);

		TR.Debug << "Total score " << total_score << std::endl;
		for ( core::Size i = 1; i <= rdc_data_all.get_number_alignment_media(); ++i  ) {
			TR.Debug << "Score for alignment medium " << rdc_data_all.get_rdc_multiset_vec()[i]->get_alignment_medium_label() << ": " << scores_all_media[i] << std::endl;
			TR.Debug << "Alignment tensor: " << std::endl;
			RDCTensorOP tensor(utility::pointer::const_pointer_cast< RDCTensor >(tensors_all_media[i]));
			if ( rdc_data_all.get_rdc_multiset_vec()[i]->get_computation_type() == RDCMultiSet::SVD ) {
				tensor->set_Dmax(21584.630);
				tensor->diagonalize_tensor();
				tensor->reorder_tensor();
			}
			tensor->show_tensor_stats(TR.Debug);
		}

		TS_ASSERT_DELTA(tensors_all_media[1]->get_ax(), -7.180E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_rh(), -6.250E-05, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_Da(), -7.749, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_R() ,  0.087, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[2]->get_ax(), -7.316E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_rh(), -1.223E-17, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_Da(), -7.760, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_R() ,  0.100, 1e-1);
	}

	/// @brief Test RDC score calculation by calling RDCData's compute_score_all_media() method.
	///        Here, the symmetric spins are explicitly declared in the input file
	///        and an identical set of rdcs is assigned the symmetric subunit.
	void test_RDCData_score_calculation_all_media_manual_symmetry_handling() {
		using namespace core::scoring::nmr::rdc;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type NH -nmr:rdc:multiset_weights 2.0 1.5");

		TR << "Testing RDCData score calculation for symmetric pose" << std::endl;

		// Setup RDC data
		RDCData rdc_data_all("core/scoring/nmr/rdc/rdc_data_2jwk_all_su_input_file.txt", tolR_);

		utility::vector1< core::Real > scores_all_media;
		utility::vector1< RDCTensorCOP > tensors_all_media;
		core::Real total_score = rdc_data_all.compute_score_all_media(tolR_, scores_all_media, tensors_all_media);

		TR.Debug << "Total score " << total_score << std::endl;
		for ( core::Size i = 1; i <= rdc_data_all.get_number_alignment_media(); ++i  ) {
			TR.Debug << "Score for alignment medium " << rdc_data_all.get_rdc_multiset_vec()[i]->get_alignment_medium_label() << ": " << scores_all_media[i] << std::endl;
			TR.Debug << "Alignment tensor: " << std::endl;
			RDCTensorOP tensor(utility::pointer::const_pointer_cast< RDCTensor >(tensors_all_media[i]));
			if ( rdc_data_all.get_rdc_multiset_vec()[i]->get_computation_type() == RDCMultiSet::SVD ) {
				tensor->set_Dmax(21584.630);
				tensor->diagonalize_tensor();
				tensor->reorder_tensor();
			}
			tensor->show_tensor_stats(TR.Debug);
		}

		TS_ASSERT_DELTA(tensors_all_media[1]->get_ax(), -7.180E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_rh(), -6.250E-05, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_Da(), -7.749, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_R() ,  0.087, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[2]->get_ax(), -7.316E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_rh(), -1.223E-17, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_Da(), -7.760, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_R() ,  0.100, 1e-1);
	}
};
