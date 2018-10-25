// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/indexed_structure_store/DirectSegmentLookup.cxxtest.hh
/// @brief Test suite for indexed_structure_store DirectSegmentLookup components.
/// @author Alex Ford (fordas@uw.edu)
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
#include <core/pose/PDBInfo.hh>
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

#include <test/protocols/indexed_structure_store/utility.hh>

#include <protocols/indexed_structure_store/DirectSegmentLookup.hh>
#include <protocols/indexed_structure_store/movers/DirectSegmentLookupMover.hh>

static basic::Tracer TR("protocols.indexed_structure_store.DirectSegmentLookup.cxxtest");

namespace protocols { namespace indexed_structure_store {
void to_json(json& j, DirectSegmentLookupResult const & r) {
	j = json{
		{"result_residues", r.result_residues},
		{"query_results", r.query_results},
		};
}
} }

namespace
{
using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::movers;
using namespace protocols::indexed_structure_store::search;

class DirectSegmentLookupTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		protocols_init();
	}

	void test_simple_mover() {
		std::string store_path = "protocols/indexed_structure_store/test_structure_store.json";
		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

		DirectSegmentLookupConfig config;
		config.rmsd_tolerance = .5;
		config.segment_cluster_tolerance = 2;
		config.max_insertion_length = 7;

		DirectSegmentLookupMover mover;
		mover.structure_store_path(store_path);
		mover.lookup_config(config);
		mover.stored_subset_name("unittest");
		mover.label_insertion("unittest_label");

		auto test_residues = store->get_residues(5);

		int tst_start = 149;
		int tst_end = 154;

		ndarray::Array<ResidueEntry, 1, 1> query_struct( test_residues.getSize<0>() - (tst_end - tst_start) );
		query_struct[ndarray::view(0,tst_start)] = test_residues[ndarray::view(0, tst_start)];
		query_struct[ndarray::view(tst_start,query_struct.getSize<0>())] = test_residues[ndarray::view(tst_end, test_residues.getSize<0>())];
		for ( core::Size i = 0; i < query_struct.getSize<0>(); ++i ) {
			query_struct[i].chain_ending = false;
		}
		query_struct[tst_start - 1].chain_ending = true;

		core::pose::PoseOP query_pose = residue_entries_to_pose(query_struct);

		TS_ASSERT_EQUALS(query_pose->num_jump(), 1);
		TS_ASSERT_EQUALS(query_pose->conformation().num_chains(), 2);
		TS_ASSERT_EQUALS(query_pose->size(), query_struct.getSize<0>());

		mover.apply(*query_pose);
		TS_ASSERT_EQUALS(query_pose->num_jump(), 0);
		TS_ASSERT_EQUALS(query_pose->conformation().num_chains(), 1);
		TS_ASSERT_EQUALS(query_pose->size(), test_residues.getSize<0>());

		core::pose::PoseOP test_pose = residue_entries_to_pose( test_residues );
		for ( core::Size i = 2; i <= test_pose->total_residue(); ++i ) {
			TS_ASSERT_EQUALS(test_pose->residue(i).has_lower_connect(), true);
		}
		for ( core::Size i = 1; i <= test_pose->total_residue() - 1; ++i ) {
			TS_ASSERT_EQUALS(test_pose->residue(i).has_upper_connect(), true);
		}

		for ( core::Size i = 2; i <= query_pose->total_residue(); ++i ) {
			TS_ASSERT_EQUALS(query_pose->residue(i).has_lower_connect(), true);
		}
		for ( core::Size i = 1; i <= query_pose->total_residue() - 1; ++i ) {
			TS_ASSERT_EQUALS(query_pose->residue(i).has_upper_connect(), true);
		}

		for ( core::Size i = 1; i <= test_pose->total_residue(); ++i ) {
			TS_ASSERT_DELTA((test_pose->residue(i).xyz("N") - query_pose->residue(i).xyz("N")).norm(), 0, 0.1);
			TS_ASSERT_DELTA((test_pose->residue(i).xyz("CA") - query_pose->residue(i).xyz("CA")).norm(), 0, 0.1);
			TS_ASSERT_DELTA((test_pose->residue(i).xyz("C") - query_pose->residue(i).xyz("C")).norm(), 0, 0.1);
			TS_ASSERT_DELTA((test_pose->residue(i).xyz("O") - query_pose->residue(i).xyz("O")).norm(), 0, 0.1);
		}

		using namespace core::select::residue_selector;
		CachedResidueSubset & cached_subsets = CachedResidueSubset::from_pose_datacache(*query_pose);
		TS_ASSERT_EQUALS(cached_subsets.has_subset("unittest"), true);
		TS_ASSERT_EQUALS(cached_subsets.get_subset("unittest")->size(), query_pose->size());
		ResidueSubset insertion_subset(*cached_subsets.get_subset("unittest"));
		for ( int r = 1; r < tst_start - 2; ++r ) {
			TS_ASSERT_EQUALS((bool)insertion_subset[r], false);
		}
		for ( int r = tst_start - 2 + 1; r < tst_end + 2 + 1; ++r ) {
			TS_ASSERT_EQUALS((bool)insertion_subset[r], true);
		}
		for ( int r = tst_end + 2 + 1; r < (int)insertion_subset.size(); ++r ) {
			TS_ASSERT_EQUALS((bool)insertion_subset[r], false);
		}

		for ( int r = 1; r <= (int)insertion_subset.size(); ++r ) {
			TS_ASSERT_EQUALS(
				(bool)query_pose->pdb_info()->res_haslabel(r, "unittest_label"),
				(bool)insertion_subset[r]
			);
		}

	}

	void test_simple_sample() {
		std::string store_path = "protocols/indexed_structure_store/test_structure_store.json";
		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

		auto test_residues = store->get_residues(5);

		int tst_start = 149;
		int tst_end = 154;
		for ( int i = tst_start - 10; i < tst_end + 10; ++i ) {
			TS_ASSERT_EQUALS( test_residues[i].chain_ending, false);
		}

		ndarray::Array<ResidueEntry, 1, 1> query_struct( test_residues.getSize<0>() - (tst_end - tst_start) );
		query_struct[ndarray::view(0,tst_start)] = test_residues[ndarray::view(0, tst_start)];
		query_struct[ndarray::view(tst_start,query_struct.getSize<0>())] = test_residues[ndarray::view(tst_end, test_residues.getSize<0>())];

		core::pose::PoseOP work_pose = residue_entries_to_pose(query_struct, "fa_standard", false);

		StructureDatabase search_db;
		search_db.initialize( store->residue_entries );

		DirectSegmentLookup::Config config;
		config.rmsd_tolerance = .5;
		config.segment_cluster_tolerance = 2;
		config.max_insertion_length = tst_end - tst_start;
		DirectSegmentLookup lookup(config);

		auto lookup_results = lookup.segment_lookup(
			store->residue_entries,
			search_db,
			*work_pose,
			tst_start - 1, tst_start + 1,
			tst_start + 1, tst_start + 3
		);

		TS_ASSERT_EQUALS(lookup_results.size(), 1);
		if ( lookup_results.size() < 1 ) {
			return;
		}

		TS_ASSERT_EQUALS(
			lookup_results.front().query_results[0].fragment_a_structure_start, tst_start - 2);
		TS_ASSERT_EQUALS(
			lookup_results.front().query_results[0].fragment_b_structure_start, tst_end);
	}

	void test_read_store()
	{
		std::string aas = "ACDEFGHIKLMNPQRSTVWY";

		std::string store_path = "protocols/indexed_structure_store/test_structure_store.json";
		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

		ndarray::Array<ResidueEntry, 1, 1> test_residues = ndarray::copy(store->get_residues(9));
		for ( auto & r: test_residues ) {
			r.structure_id = 0;
			TS_ASSERT_DIFFERS( aas.find(r.sc.aa), std::string::npos );
		}
		test_residues[test_residues.getSize<0>() - 1].chain_ending = true;

		core::pose::PoseOP work_pose = residue_entries_to_pose(
			store->get_residues(9), "fa_standard", true);

		TS_ASSERT_EQUALS(test_residues.getSize<0>(), work_pose->total_residue());

		auto extracted_residues = extract_residue_entries(*work_pose, false);
		assert_residue_entries_almost_equal(test_residues, extracted_residues);

		core::pose::PoseOP work_pose2 = residue_entries_to_pose(
			extracted_residues, "fa_standard", true);
		auto extracted_residues2 = extract_residue_entries(*work_pose2, false);

		assert_residue_entries_almost_equal( extracted_residues, extracted_residues2 );
	}
};

}
