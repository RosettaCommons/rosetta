// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREEnergy.cxxtest.hh
/// @brief   unit test for class PREEnergy which is the energy method to score with NMR PRE data
/// @details Last modified: 12/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>

// Unit headers
#include <protocols/nmr/pre/PREEnergy.hh>
#include <core/scoring/nmr/pre/PREData.hh>
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

static basic::Tracer TR("protocols.nmr.pre.PREEnergy.cxxtest");

class PREEnergyTests : public CxxTest::TestSuite {

private:
	core::pose::Pose gb1_;

public:

	/// @brief Setup Test
	void setUp() {
		// Init options from pseudo command line
		core_init();

		// Load pose from pdb
		core::import_pose::pose_from_file(gb1_, "protocols/nmr/pre/gb1_.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		gb1_.clear();
	}

	void test_PREEnergy_instantiation() {
		using namespace core::scoring::nmr::pre;
		using namespace core::scoring;
		using namespace protocols::nmr::pre;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:input_file protocols/nmr/pre/pre_data_input.txt -nmr:pre:multiset_weights 1.0 -nmr:pre:show_info true");

		TR << "Test instantiation of PREEnergy." << std::endl;

		PREEnergy pre_energy;
		PREData & pre_data = pre_energy.get_pre_data_from_pose(gb1_);
		pre_data.show(TR.Debug);
		core::Size number_spinlabel_sites = pre_data.get_number_spinlabel_sites();
		TS_ASSERT_EQUALS(number_spinlabel_sites, 1);
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pre , 2.0 );
		TS_ASSERT(sfxn.has_nonzero_weight( nmr_pre ));
	}

	void test_score_calculation_with_spinlabel_distance_filter() {
		using namespace core::scoring::nmr::pre;
		using namespace core::scoring;
		using namespace protocols::nmr::pre;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true -nmr:pre:input_file protocols/nmr/pre/pre_data_input.txt -nmr:pre:multiset_weights 1.0 -nmr:pre:show_info true -nmr:spinlabel:highres_conformer_filter_type DISTANCE");

		TR << "Test PREEnergy calculation with implicit spinlabel." << std::endl;

		PREEnergy pre_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pre , 2.0 );
		EnergyMap emap;
		pre_energy.finalize_total_energy(gb1_, sfxn, emap);
		TR.Debug << "PRE score calculated with spinlabel ensemble and neighbor count filter: " << emap[ nmr_pre ] << std::endl;
		TS_ASSERT_DELTA(emap[ nmr_pre ], 16.4859, 1e-1 );
	}

	void test_score_calculation_with_spinlabel_bump_energy_filter() {
		using namespace core::scoring::nmr::pre;
		using namespace core::scoring;
		using namespace protocols::nmr::pre;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true -nmr:pre:input_file protocols/nmr/pre/pre_data_input.txt -nmr:pre:multiset_weights 1.0 -nmr:pre:show_info true -nmr:spinlabel:highres_conformer_filter_type BUMP_ENERGY");

		TR << "Test PREEnergy calculation with explicit spinlabel and no side chain packing." << std::endl;

		PREEnergy pre_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pre , 2.0 );
		EnergyMap emap;
		pre_energy.finalize_total_energy(gb1_, sfxn, emap);
		TR.Debug << "PRE score calculated with dummy spinlabel ensemble and bump energy filter: " << emap[ nmr_pre ] << std::endl;
		TS_ASSERT_DELTA(emap[ nmr_pre ], 24.7708, 1e-1 );
	}

	void test_show_info() {
		using namespace core::scoring::nmr::pre;
		using namespace core::scoring;
		using namespace protocols::nmr::pre;
		using namespace core::pose::datacache;

		// Additional options; sets also fixed RG seed for NLS fitting
		core_init_with_additional_options("-nmr:pre:normalize_data true -nmr:pre:input_file protocols/nmr/pre/pre_data_input.txt -nmr:pre:multiset_weights 1.0 -nmr:pre:show_info true -nmr:spinlabel:highres_conformer_filter_type BUMP_ENERGY");

		TR << "Test PREEnergy show info function." << std::endl;

		PREEnergy pre_energy;
		ScoreFunction sfxn;
		sfxn.set_weight( nmr_pre , 2.0 );
		EnergyMap emap;
		pre_energy.finalize_total_energy(gb1_, sfxn, emap);
		pre_energy.show_additional_info(TR.Debug, gb1_, true);
		// The function above creates a PREdata object and puts it into the pose's datacache
		TS_ASSERT(gb1_.data().has( CacheableDataType::NMR_PRE_DATA ));
	}
};
