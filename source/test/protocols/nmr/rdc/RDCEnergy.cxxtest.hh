// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/rdc/RDCEnergy.cxxtest.hh
/// @brief   unit test for class RDCEnergy which is the energy method to score the pose with NMR RDC data
/// @details Last modified: 12/09/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

// Unit headers
#include <protocols/nmr/rdc/RDCEnergy.hh>
#include <core/scoring/nmr/rdc/RDCData.hh>
#include <core/scoring/nmr/NMRDataFactory.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/datacache/CacheableDataType.hh>

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

static basic::Tracer TR("protocols.nmr.rdc.RDCEnergy.cxxtest");

class RDCEnergyTests : public CxxTest::TestSuite {

private: // Data
	core::pose::Pose g_;

public: // Methods

	/// @brief Setup Test
	void setUp() {
		// Init options from pseudo command line
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type none -nmr:rdc:input_file protocols/nmr/rdc/rdc_data_g_input_file.txt -nmr:rdc:multiset_weights 1.0 1.0 1.0 1.0 -nmr:rdc:show_info true");
		// Load pose from pdb
		core::import_pose::pose_from_file(g_, "protocols/nmr/rdc/g_nmr_A.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		g_.clear();
	}

	void test_RDCEnergy_instantiation() {
		using namespace core::scoring::nmr::rdc;
		using namespace core::scoring;
		using namespace protocols::nmr::rdc;

		RDCEnergy rdc_energy;
		RDCData & rdc_data = rdc_energy.get_rdc_data_from_pose(g_);
		rdc_data.show(TR.Debug);
		core::Size number_alignment_media = rdc_data.get_number_alignment_media();
		TS_ASSERT_EQUALS(number_alignment_media, 2);
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_rdc , 2.0 );
		TS_ASSERT(sfxn.has_nonzero_weight( nmr_rdc ));
	}

	void test_RDCEnergy_score_calculation() {
		using namespace core::scoring::nmr::rdc;
		using namespace core::scoring;
		using namespace protocols::nmr::rdc;

		// Set fixed RG seed for NLS fitting
		initialize_rng();

		RDCEnergy rdc_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_rdc , 2.0 );
		EnergyMap emap;
		rdc_energy.finalize_total_energy(g_, sfxn, emap);
		TS_ASSERT_DELTA(emap[ nmr_rdc ], 1564.388, 1e-1 );
	}

	void test_PCSEnergy_show_info() {
		using namespace protocols::nmr::rdc;
		using namespace core::pose::datacache;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:rdc:correct_sign false -nmr:rdc:normalization_type none -nmr:rdc:input_file protocols/nmr/rdc/rdc_data_g_input_file.txt -nmr:rdc:multiset_weights 1.0 1.0 1.0 1.0 -nmr:rdc:show_info true -out:level 500");

		RDCEnergy rdc_energy;
		rdc_energy.show_additional_info(TR.Debug, g_); // performs score calculation and writes out tensor params
		// The function above creates a RDCdata object and puts it into the pose's datacache
		TS_ASSERT(g_.data().has( CacheableDataType::NMR_RDC_DATA ));
	}
};
