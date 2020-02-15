// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/HbondMetricTests.cxxtest.hh
/// @brief  Unit test for hbond metric.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/per_residue_metrics/HbondMetric.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/util.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("HbondMetricTests");

using namespace core::simple_metrics;
using namespace core::simple_metrics::per_residue_metrics;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace core::select;
using namespace core::scoring::hbonds;

class HbondMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		//PDB sugars, ignore unrecognized, don't load pdb components.
		// Skip waters. Use beta scorefunction for hbonding corrections.
		std::string opts =
			"-include_sugars "
			"-auto_detect_glycan_connections "
			"-alternate_3_letter_codes pdb_sugar "
			"-maintain_links "
			"-ignore_waters true "
			"-beta_nov16 "
			"-corrections::gen_potential "
			"-include_vrt false "
			"-ignore_unrecognized_res "
			"-ignore_zero_occupancy false "
			"-load_PDB_components false ";

		core_init_with_additional_options(opts);
		core::import_pose::pose_from_file( pose_, "core/simple_metrics/2cl2_opt_waters2.pdb", false, core::import_pose::PDB_file );

		selector_ = utility::pointer::make_shared< GlycanResidueSelector >();
		not_selector_ = utility::pointer::make_shared< NotResidueSelector >(selector_);
		true_selector_ = utility::pointer::make_shared< TrueResidueSelector >();

	}

	void tearDown(){

	}

	void run_test( bool include_self){
		std::map< core::Size, core::Real > hbonds;
		HBondSet hb_set = HBondSet(pose_, false);

		utility::vector1< core::Size > sele1 = get_residues_from_subset(selector_->apply(pose_));
		utility::vector1< core::Size > sele2 = get_residues_from_subset(not_selector_->apply(pose_));

		HbondMetric metric = HbondMetric();
		metric.set_residue_selector(selector_);
		metric.set_residue_selector2(not_selector_);
		metric.set_include_self(include_self);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		for ( core::Size const s1 : sele1 ) {
			utility::vector1< HBondCOP > s1_hbonds = hb_set.residue_hbonds(s1);
			for ( HBondCOP hb : s1_hbonds ) {
				core::Size s2 = next_hb_res(*hb, s1);
				if ( sele2.contains(s2) ) {
					if ( s2 == s1 && ! include_self ) {
						continue;
					}

					TS_ASSERT(per_res.count(s1));
					hbonds[s1]+=1;
				}
			}
		}

		for ( auto const & p : hbonds ) {
			TS_ASSERT_EQUALS(per_res.at(p.first), p.second);
		}
	}
	void test_default(){
		run_test(false);
	}

	void test_include_self(){
		run_test(true);
	}


private:

	core::pose::Pose pose_;
	core::select::residue_selector::GlycanResidueSelectorCOP selector_;
	core::select::residue_selector::NotResidueSelectorCOP not_selector_;
	core::select::residue_selector::TrueResidueSelectorCOP true_selector_;




};
