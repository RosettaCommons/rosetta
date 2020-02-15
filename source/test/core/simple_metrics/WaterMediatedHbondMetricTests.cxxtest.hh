// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/simple_metrics/WaterMediatedHbondMetricTests.cxxtest.hh
/// @brief  Unit test for water mediated hbond metric.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/util.hh>

#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("WaterMediatedHbondMetricTests");

using namespace core::simple_metrics;
using namespace core::simple_metrics::per_residue_metrics;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace core::select;
using namespace core::scoring::hbonds;

class WaterMediatedHbondMetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){

		//PDB sugars, ignore unrecognized, don't load pdb components.
		// Load waters. Use beta scorefunction for hbonding corrections.
		std::string opts =
			"-include_sugars "
			"-auto_detect_glycan_connections "
			"-alternate_3_letter_codes pdb_sugar "
			"-maintain_links "
			"-ignore_waters false "
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
		index_selector_ = utility::pointer::make_shared<  ResidueIndexSelector >();

	}

	void tearDown(){

	}

	void test_basic_functions(){

		//Test our main recursive function that identifies paths for each residue
		// Paths came from function prototyped and tested first in PyRosetta
		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();
		HBondSet hb_set = HBondSet(pose_, false);

		// Populate all waters
		utility::vector1< core::Size > all_waters;
		for ( core::Size i = 1; i <= pose_.size(); ++i ) {
			if ( pose_.residue_type(i).is_water() ) {
				all_waters.push_back(i);
			}
		}

		//Depth one
		std::map< std::string, core::Size > paths_depth_one;
		utility::vector1< std::string > paths_304;
		paths_304.push_back("304-467-50");
		paths_304.push_back("304-819-290");

		metric.find_hb_paths(hb_set, all_waters, paths_depth_one, 304, 1);

		TS_ASSERT(paths_depth_one.size() == paths_304.size());
		for ( std::string const & s : paths_304 ) {
			TS_ASSERT(paths_depth_one.count(s));
		}

		//Depth 1 + 2
		std::map< std::string, core::Size > paths_depth_two;
		paths_304.push_back("304-467-475-51");
		paths_304.push_back("304-467-475-302");
		paths_304.push_back("304-467-833-302");
		paths_304.push_back("304-834-688-183");
		paths_304.push_back("304-834-818-290");
		paths_304.push_back("304-837-566-90");

		metric.find_hb_paths(hb_set, all_waters, paths_depth_two, 304, 2);
		TS_ASSERT(paths_depth_two.size() == paths_304.size());
		for ( std::string const & s : paths_304 ) {
			TS_ASSERT(paths_depth_two.count(s));
		}

		//Depth 1 + 2 + 3
		std::map< std::string, core::Size > paths_depth_three;
		paths_304.push_back("304-467-833-845-303");
		paths_304.push_back("304-834-424-567-90");
		paths_304.push_back("304-834-818-817-288");

		metric.find_hb_paths(hb_set, all_waters, paths_depth_three, 304, 3);
		TS_ASSERT(paths_depth_three.size() == paths_304.size());
		for ( std::string const & s : paths_304 ) {
			TS_ASSERT(paths_depth_three.count(s));
		}

		// Assert that all bridge residues are waters
		for ( auto const & p : paths_depth_three ) {
			utility::vector1< std::string > SP = utility::string_split(p.first, '-');
			for ( core::Size i=1; i <= (SP.size()-1); ++i ) {
				if ( i == 1 ) continue;

				TS_ASSERT(pose_.residue_type(utility::string2Size(SP[i])).is_water());
			}
		}
	}

	void test_depth_one() {
		//Test our main recursive function that identifies paths for each residue
		// Paths came from function prototyped and tested first in PyRosetta
		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();
		metric.set_residue_selector(selector_);

		//Use the true selector to get ALL bridges between glycans and Everything else, including self and other glycan residues.
		metric.set_residue_selector2(true_selector_);

		//These are default - but just to make sure they are set.
		metric.set_include_self(false);
		metric.set_depth(1);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		//Assert we only have the residues in the selector here.
		TS_ASSERT(per_res.size() == get_residues_from_subset(selector_->apply(pose_)).size());

		//Assert we bound the correct number of unique paths
		TS_ASSERT(per_res.at(299) == 2);
		TS_ASSERT(per_res.at(300) == 1);
		TS_ASSERT(per_res.at(301) == 1);
		TS_ASSERT(per_res.at(302) == 1);
		TS_ASSERT(per_res.at(303) == 0);
		TS_ASSERT(per_res.at(304) == 2);

		//305 makes a bridge to itself through a water
		// Without the option, we don't count it.
		TS_ASSERT(per_res.at(305) == 0);

	}

	void test_depth_two() {
		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();
		metric.set_residue_selector(selector_);

		//Use the true selector to get ALL bridges between glycans and Everything else, including self and other glycan residues.
		metric.set_residue_selector2(true_selector_);

		metric.set_include_self(false);
		metric.set_depth(2);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		//Assert we only have the residues in the selector here.
		TS_ASSERT(per_res.size() == get_residues_from_subset(selector_->apply(pose_)).size());

		TS_ASSERT(per_res.at(299) == 3);
		TS_ASSERT(per_res.at(300) == 1);
		TS_ASSERT(per_res.at(301) == 1);
		TS_ASSERT(per_res.at(302) == 6);
		TS_ASSERT(per_res.at(303) == 1);
		TS_ASSERT(per_res.at(304) == 8);
		TS_ASSERT(per_res.at(305) == 0);


	}

	void test_only_depth_two() {
		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();

		//Speed
		index_selector_->append_index(299);
		index_selector_->append_index(300);
		index_selector_->append_index(302);

		metric.set_residue_selector(index_selector_);

		//Use the true selector to get ALL bridges between glycans and Everything else, including self and other glycan residues.
		metric.set_residue_selector2(true_selector_);

		metric.set_include_self(false);
		metric.set_depth(2);
		metric.set_include_only_set_depth(true);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		TS_ASSERT(per_res.size() == get_residues_from_subset(index_selector_->apply(pose_)).size());

		TS_ASSERT(per_res.at(299) == 1);
		TS_ASSERT(per_res.at(300) == 0);
		TS_ASSERT(per_res.at(302) == 5);

	}

	void test_depth_one_with_self() {
		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();
		index_selector_->set_index(305);
		metric.set_residue_selector(index_selector_);

		//Use the true selector to get ALL bridges between glycans and Everything else, including self and other glycan residues.
		metric.set_residue_selector2(true_selector_);
		metric.set_include_self(true);
		metric.set_depth(1);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		//305 Makes a bridge to itself through a water
		// Make sure it is there as we add the option to include self.
		TS_ASSERT(per_res.at(305) == 1);
	}

	void test_sele1_to_sele2() {
		//Assert that we are only getting paths to second selector.

		WaterMediatedHbondMetric metric = WaterMediatedHbondMetric();

		//Speed
		index_selector_->append_index(299);
		index_selector_->append_index(300);
		index_selector_->append_index(302);

		metric.set_residue_selector(index_selector_);

		//Use the true selector to get ALL bridges between glycans and Everything else, including self and other glycan residues.
		metric.set_residue_selector2(not_selector_);
		metric.set_include_self(true);
		metric.set_depth(2);

		std::map< core::Size, core::Real > const per_res = metric.calculate(pose_);

		TS_ASSERT(per_res.at(299) == 3);
		TS_ASSERT(per_res.at(300) == 0);
		TS_ASSERT(per_res.at(302) == 3);

	}

private:

	core::pose::Pose pose_;
	core::select::residue_selector::GlycanResidueSelectorCOP selector_;
	core::select::residue_selector::NotResidueSelectorCOP not_selector_;
	core::select::residue_selector::TrueResidueSelectorCOP true_selector_;
	core::select::residue_selector::ResidueIndexSelectorOP index_selector_;

};
