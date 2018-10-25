// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/indexed_structure_store/SegmentSequenceProfile.cxxtest.hh
/// @brief Test suite for indexed_structure_store SegmentSequenceProfile components.
/// @author Alex Ford (fordas@uw.edu)
//
//

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <util/pose_funcs.hh>
#include <boost/range/combine.hpp>

// Project headers
#include <basic/Tracer.hh>

#include <boost/range.hpp>
#include <boost/range/adaptors.hpp>

#include <numeric/xyz.json.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/select/residue_selector/CachedResidueSubset.hh>

#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.json.hh>
#include <protocols/indexed_structure_store/utility.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.json.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>
#include <protocols/indexed_structure_store/SegmentSequenceProfile.hh>
#include <protocols/indexed_structure_store/movers/SegmentSequenceProfileMover.hh>

#include <test/protocols/indexed_structure_store/utility.hh>

static basic::Tracer TR("protocols.indexed_structure_store.SegmentSequenceProfile.cxxtest");

namespace
{
using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::movers;
using namespace protocols::indexed_structure_store::search;

class SegmentSequenceProfileTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		protocols_init();
	}

	void test_simple_profile() {
		std::string store_path = "protocols/indexed_structure_store/test_store/structure/pdb_dir";

		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

		StructureDatabaseOP database(new StructureDatabase());
		database->initialize(store->residue_entries);

		auto test_residues = store->get_residues(0);
		core::pose::PoseOP test_pose = residue_entries_to_pose(test_residues);

		SegmentSequenceProfile profiler;
		profiler.config.pseudocount = 0;

		SegmentSequenceProfileResult profile_result = profiler.segment_profile(
			*store, *database, *test_pose, 1, test_pose->total_residue() + 1);

		// Test database has four 19-residue structures, with varying last residue
		TS_ASSERT_EQUALS(profile_result.counts.rows(), 19);
		for ( core::Size i = 0; i < 18; ++i ) {
			TS_ASSERT_EQUALS(profile_result.counts.row(i).sum(), 4);
			TS_ASSERT_EQUALS(profile_result.counts(i, core::chemical::aa_from_oneletter_code(test_residues[i].sc.aa) - core::chemical::first_l_aa), 4);
		}

		TS_ASSERT_EQUALS(profile_result.counts.row(18).sum(), 4);
		TS_ASSERT_EQUALS(profile_result.counts(18, core::chemical::aa_from_oneletter_code('A') - core::chemical::first_l_aa), 1);
		TS_ASSERT_EQUALS(profile_result.counts(18, core::chemical::aa_from_oneletter_code('N') - core::chemical::first_l_aa), 1);
		TS_ASSERT_EQUALS(profile_result.counts(18, core::chemical::aa_from_oneletter_code('S') - core::chemical::first_l_aa), 1);
		TS_ASSERT_EQUALS(profile_result.counts(18, core::chemical::aa_from_oneletter_code('V') - core::chemical::first_l_aa), 1);
	}

};

}
