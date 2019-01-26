// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSEnergy.cxxtest.hh
/// @brief   unit test for class PCSEnergy which is the energy method for working with PCS data
/// @details Last modified: 07/26/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

// Unit headers
#include <protocols/nmr/pcs/PCSEnergy.hh>
#include <core/scoring/nmr/pcs/PCSData.hh>

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

static basic::Tracer TR("protocols.nmr.pcs.PCSEnergy.cxxtest");

class PCSEnergyTests : public CxxTest::TestSuite {

private:
	core::pose::Pose psR2_;

public:

	/// @brief Setup Test
	void setUp() {
		// Init options from pseudo command line
		core_init_with_additional_options("-nmr:pcs:input_file protocols/nmr/pcs/pcs_data_2ksy_input_file.txt -nmr:pcs:optimize_tensor true -nmr:pcs:nls_repeats 10 -nmr:pcs:multiset_weights 1.0 1.0 -nmr:pcs:show_info true");
		// Load pose from pdb
		core::import_pose::pose_from_file(psR2_, "protocols/nmr/pcs/2ksy_chainA_renum.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		psR2_.clear();
	}

	// AMW: Right now these tests only run on mac.clang because they are too
	// slow with STL-debug. For them to run on all four platforms, please let
	// me know and we can work through some solutions. I bet we don't need
	// four data sets from each of two spinlabel positions on a 200 residue
	// protein to get an idea of it this code is working.

	void test_PCSEnergy_instantiation() {
#ifdef _GLIBCXX_DEBUG
		TS_ASSERT( true );
#else
		using namespace core::scoring::nmr::pcs;
		using namespace core::scoring;
		using namespace protocols::nmr::pcs;

		PCSEnergy pcs_energy;
		PCSData & pcs_data = pcs_energy.get_pcs_data_from_pose(psR2_);
		pcs_data.show(TR.Debug);
		core::Size number_tags = pcs_data.get_number_tags();
		TS_ASSERT_EQUALS(number_tags, 2);
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pcs , 2.0 );
		TS_ASSERT(sfxn.has_nonzero_weight( nmr_pcs ));
#endif
	}

	void test_PCSEnergy_score_calculation() {
#ifdef _GLIBCXX_DEBUG
		TS_ASSERT( true );
#else
		using namespace core::scoring::nmr::pcs;
		using namespace core::scoring;
		using namespace protocols::nmr::pcs;

		// Set fixed RG seed for NLS fitting
		initialize_rng();

		PCSEnergy pcs_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pcs , 2.0 );
		EnergyMap emap;
		pcs_energy.finalize_total_energy(psR2_, sfxn, emap);
		TS_ASSERT_DELTA(emap[ nmr_pcs ], 0.000635452, 1e-3 );
#endif
	}

	void test_PCSEnergy_show_info() {
#ifdef _GLIBCXX_DEBUG
		TS_ASSERT( true );
#else
		using namespace protocols::nmr::pcs;
		using namespace core::pose::datacache;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pcs:input_file protocols/nmr/pcs/pcs_data_2ksy_input_file.txt -nmr:pcs:optimize_tensor true -nmr:pcs:nls_repeats 10 -nmr:pcs:multiset_weights 1.0 1.0 -nmr:pcs:show_info true -out:level 500");

		PCSEnergy pcs_energy;
		pcs_energy.show_additional_info(TR.Debug, psR2_); // performs score calculation and writes out tensor params
		// The function above creates a PCSdata object and puts it into the pose's datacache
		TS_ASSERT(psR2_.data().has( CacheableDataType::NMR_PCS_DATA ));
#endif
	}

	void test_PCSEnergy_score_calculation_scaled_data() {
#ifdef _GLIBCXX_DEBUG
		TS_ASSERT( true );
#else
		using namespace core::scoring::nmr::pcs;
		using namespace core::scoring;
		using namespace protocols::nmr::pcs;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pcs:input_file protocols/nmr/pcs/pcs_data_2ksy_input_file.txt -nmr:pcs:normalize_data true -nmr:pcs:optimize_tensor true -nmr:pcs:nls_repeats 10 -nmr:pcs:multiset_weights 1.0 1.0");

		PCSEnergy pcs_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pcs , 2.0 );
		EnergyMap emap;
		pcs_energy.finalize_total_energy(psR2_, sfxn, emap);
		TS_ASSERT_DELTA(emap[ nmr_pcs ], 0.000966956, 1e-3 );
#endif
	}

	void test_PCSEnergy_show_info_scaled_data() {
#ifdef _GLIBCXX_DEBUG
		TS_ASSERT( true );
#else
		using namespace protocols::nmr::pcs;
		using namespace core::pose::datacache;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pcs:input_file protocols/nmr/pcs/pcs_data_2ksy_input_file.txt -nmr:pcs:normalize_data true -nmr:pcs:optimize_tensor true -nmr:pcs:nls_repeats 10 -nmr:pcs:multiset_weights 1.0 1.0 -nmr:pcs:show_info true -out:level 500");

		PCSEnergy pcs_energy;
		pcs_energy.show_additional_info(TR.Debug, psR2_); // performs score calculation and writes out tensor params
		// The function above creates a PCSdata object and puts it into the pose's datacache
		TS_ASSERT(psR2_.data().has( CacheableDataType::NMR_PCS_DATA ));
#endif
	}

};
