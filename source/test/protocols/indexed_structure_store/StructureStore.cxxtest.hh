// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/indexed_structure_store/StructureStore.cxxtest.hh
/// @brief Test suite for protocols/indexed_structure_store/*
/// @author Alex Ford (fordas@uw.edu)

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>

#include <boost/range/combine.hpp>

#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.hh>
#include <protocols/indexed_structure_store/StructureStoreProvider.hh>
#include <protocols/indexed_structure_store/JSONStructureStoreBackend.hh>
#include <protocols/indexed_structure_store/utility.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.json.hh>

#include <test/protocols/indexed_structure_store/utility.hh>

static basic::Tracer TR("core.indexed_structure_store.StructureStore.cxxtest");

namespace
{
using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::search;

class StructureStoreTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		core_init();
	}

	void test_read_h5_store()
	{
#ifdef USEHDF5
		std::cout << "Load store." << std::endl;
		std::string store_path = "protocols/indexed_structure_store/test_structure_store.json";
		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

		std::cout << "Initialize search db." << std::endl;
		StructureDatabase structure_db;
		structure_db.initialize(store->residue_entries);

		std::cout << "Get residues." << std::endl;
		auto test_residues = store->get_residues(1);
		TR << std::vector<ResidueEntry>(test_residues.begin(), test_residues.end()).size() << std::endl;
#endif
	}

	void test_convert_stores() {
		std::string input_store = "protocols/indexed_structure_store/test_store/structure/pdb_dir";

		TR << "Load store: " << input_store << std::endl;
		StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store(input_store);

		std::string output_store = "protocols/indexed_structure_store/test_store/structure/test.json";
		StructureStoreManager::get_instance()->write_structure_store(output_store, *store);

		StructureStoreOP restore = StructureStoreManager::get_instance()->load_structure_store(output_store);
		assert_residue_entries_almost_equal(store->residue_entries, restore->residue_entries);
	}

	void test_read_stores()
	{
		std::vector<std::string> test_stores =  {
			"protocols/indexed_structure_store/test_store/structure/pdb_dir",
			"protocols/indexed_structure_store/test_store/structure/silent_store.silent",
			};

		for ( std::string store_path : test_stores ) {

			TR << "Load store: " << store_path << std::endl;
			StructureStoreOP store = StructureStoreManager::get_instance()->load_structure_store( store_path );

			TR.Debug << "Initialize search db." << std::endl;
			StructureDatabase structure_db;
			structure_db.initialize(store->residue_entries);

			TR.Debug << "Get structures." << std::endl;
			TR.Debug << json(std::vector<StructureEntry>(
				store->structure_entries.begin(), store->structure_entries.end())) << std::endl;

			TR.Debug << "Get residues." << std::endl;
			TS_ASSERT_EQUALS(store->structure_entries.getSize<0>(), 4);
			std::vector<ndarray::Array<const ResidueEntry, 1>> result_residues;
			for ( core::Size i = 0; i < 4; ++i ) {
				result_residues.push_back(store->get_residues(i));
				TS_ASSERT_EQUALS(result_residues[i].getSize<0>(), 19);
			}

			assert_residue_entries_almost_equal(
				result_residues[0][ndarray::view(0, 18)], result_residues[1][ndarray::view(0, 18)], false);
			assert_residue_entries_almost_equal(
				result_residues[0][ndarray::view(0, 18)], result_residues[2][ndarray::view(0, 18)], false);
			assert_residue_entries_almost_equal(
				result_residues[0][ndarray::view(0, 18)], result_residues[3][ndarray::view(0, 18)], false);

			assert_residue_entries_almost_equal(
				result_residues[0], result_residues[1], false, true, false);
			assert_residue_entries_almost_equal(
				result_residues[0], result_residues[2], false, true, false);
			assert_residue_entries_almost_equal(
				result_residues[0], result_residues[3], false, true, false);
		}
	}

	class MockStore : public StructureStoreProvider
	{
	public:
		virtual bool can_load(std::string path) {
			return path.substr(0, 4) == "mock";
		}

		virtual bool can_write(std::string path) {
			return path.substr(0, 4) == "mock";
		}

		virtual StructureStoreOP load_store(std::string path) {
			load_calls[path] += 1;
			return StructureStoreOP(new StructureStore());
		}

		virtual void write_store(std::string path, StructureStore &) {
			write_calls[path] += 1;
		}

		std::map<std::string, int> load_calls;
		std::map<std::string, int> write_calls;
	};

	void test_store_cache()
	{
		// Store caches databases, returning references to existing database on
		// repeated loads for a target db.
		std::shared_ptr<MockStore> mock(new MockStore);
		StructureStoreManager::get_instance()->register_store_provider(
			-1, "mock", StructureStoreProviderOP(mock));

		StructureStoreOP one_a =
			StructureStoreManager::get_instance()->load_structure_store("mock1");
		StructureStoreOP one_b =
			StructureStoreManager::get_instance()->load_structure_store("mock1");
		StructureStoreOP two_a =
			StructureStoreManager::get_instance()->load_structure_store("mock2");

		TS_ASSERT(one_a == one_b);
		TS_ASSERT(one_a != two_a);
		TS_ASSERT(mock->load_calls["mock1"] == 1);
		TS_ASSERT(mock->load_calls["mock2"] == 1);

		one_a.reset();
		one_b.reset();
		StructureStoreOP one_c =
			StructureStoreManager::get_instance()->load_structure_store("mock1");
		StructureStoreOP two_b =
			StructureStoreManager::get_instance()->load_structure_store("mock2");
		TS_ASSERT(two_a == two_b);
		TS_ASSERT(mock->load_calls["mock1"] == 1);
		TS_ASSERT(mock->load_calls["mock2"] == 1);
	}
};

}
