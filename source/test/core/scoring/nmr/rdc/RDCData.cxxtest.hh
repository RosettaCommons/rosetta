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

static basic::Tracer TR("core.scoring.nmr.rdc.RDCData.cxxtest");

class RDCDataTests : public CxxTest::TestSuite {

private: // Data
	core::pose::Pose g_;

public: // Methods

	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		// RDC data are not pre-scaled yet and we have to correct the sign of the 15N gyromagnetic ratio
		core_init_with_additional_options("-nmr:rdc:correct_sign true -nmr:rdc:normalization_type NH -nmr:rdc:multiset_weights 2.0 1.0 1.5 1.0 2.5 0.5 1.0");

		// Load pose from pdb
		core::import_pose::pose_from_file(g_, "core/scoring/nmr/rdc/gb1.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		g_.clear();
	}

	/// Test the creation of RDCData
	void test_RDCData_instantiation() {
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCData instantiation" << std::endl;

		RDCData rdc_data_all("core/scoring/nmr/rdc/rdc_data_g_input_file.txt", g_);
		rdc_data_all.show(TR.Debug);
		TS_ASSERT_EQUALS(rdc_data_all.get_number_alignment_media(), 7);
	}

	/// Test the creation of RDCData through the common interface of the NMRDataFactory.
	void test_RDCData_creation_by_NMRDataFactory() {
		using namespace core::scoring::nmr;
		using namespace core::scoring::nmr::rdc;

		TR << "Testing RDCData instantiation" << std::endl;

		RDCDataOP rdc_data_all_ptr(utility::pointer::static_pointer_cast< RDCData >( NMRDataFactory::get_instance()->get_nmr_data("RDC", "core/scoring/nmr/rdc/rdc_data_g_input_file.txt", g_) ) );
		rdc_data_all_ptr->show(TR.Debug);
		TS_ASSERT_EQUALS(rdc_data_all_ptr->get_number_alignment_media(), 7);
	}

	/// Test RDC score calculation by calling RDCData's compute_score_all_media() method
	void test_RDCData_score_calculation_all_media() {
		using namespace core::scoring::nmr::rdc;

		// Set fixed RG seed for NLS fitting
		initialize_rng();

		TR << "Testing RDCData score calculation" << std::endl;

		RDCData rdc_data_all("core/scoring/nmr/rdc/rdc_data_g_input_file.txt", g_);

		utility::vector1< core::Real > scores_all_media;
		utility::vector1< RDCTensorCOP > tensors_all_media;
		core::Real total_score = rdc_data_all.compute_score_all_media(g_, scores_all_media, tensors_all_media);

		TR.Debug << "Total score " << total_score << std::endl;
		for ( core::Size i = 1; i <= rdc_data_all.get_number_alignment_media(); ++i  ) {
			TR.Debug << "Score for alignment medium " << rdc_data_all.get_rdc_multiset_vec()[i]->get_alignment_medium_label() << ": " << scores_all_media[i] << std::endl;
			TR.Debug << "Alignment tensor: " << std::endl;
			RDCTensorOP tensor(utility::pointer::const_pointer_cast< RDCTensor >(tensors_all_media[i]));
			if ( rdc_data_all.get_rdc_multiset_vec()[i]->get_computation_type() == RDCMultiSet::SVD && !rdc_data_all.get_rdc_multiset_vec()[i]->tensor_fixed() ) {
				tensor->set_Dmax(21584.630);
				tensor->diagonalize_tensor();
				tensor->reorder_tensor();
			}
			tensor->show_tensor_stats(TR.Debug);
		}

		TS_ASSERT_DELTA(tensors_all_media[1]->get_ax(), -7.927E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_rh(), -2.318E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_Da(), -8.555, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[1]->get_R() ,  0.292, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[2]->get_ax(), -5.490E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_rh(), -3.341E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_Da(), -6.123, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[2]->get_R() ,  0.597, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[3]->get_ax(), -7.927E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[3]->get_rh(), -2.318E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[3]->get_Da(), -8.555, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[3]->get_R() ,  0.292, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[4]->get_ax(), -5.406E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[4]->get_rh(), -3.590E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[4]->get_Da(), -6.022, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[4]->get_R() ,  0.664, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[5]->get_ax(), -7.927E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[5]->get_rh(), -2.318E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[5]->get_Da(), -8.555, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[5]->get_R() ,  0.292, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[6]->get_ax(),  5.122E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[6]->get_rh(),  3.401E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[6]->get_Da(),  5.528, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[6]->get_R() ,  0.664, 1e-1);

		TS_ASSERT_DELTA(tensors_all_media[7]->get_ax(), -7.927E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[7]->get_rh(), -2.318E-04, 1e-4);
		TS_ASSERT_DELTA(tensors_all_media[7]->get_Da(), -8.555, 1e-1);
		TS_ASSERT_DELTA(tensors_all_media[7]->get_R() ,  0.292, 1e-1);
	}
};
