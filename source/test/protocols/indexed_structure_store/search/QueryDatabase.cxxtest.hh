// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/indexed_structure_store/search/QueryDatabase.cxxtest.hh
/// @brief Test suite for protocols/indexed_structure_store/search/QueryDatabase
/// @author Alex Ford (fordas@uw.edu)
//
// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

#include <vector>

// Project headers
#include <basic/Tracer.hh>

#include <protocols/indexed_structure_store/search/QueryDatabase.hh>
#include <protocols/indexed_structure_store/search/QueryDatabase.json.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.json.hh>
#include <protocols/indexed_structure_store/orient_array.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <json.hpp>

static basic::Tracer TR("protocols.indexed_structure_store.search.QueryDatabase.cxxtest");

using nlohmann::json;

namespace
{

using namespace protocols::indexed_structure_store;
using namespace protocols::indexed_structure_store::search;

class QueryDatabaseTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		core_init_with_additional_options( "" );
	}

	void test_single_query()
	{

		core::pose::PoseOP test_pose = core::import_pose::pose_from_file("protocols/indexed_structure_store/test_structure.pdb", core::import_pose::PDB_file);

		ndarray::Array<ResidueEntry, 1> test_residues = extract_residue_entries( *test_pose );

		StructureDatabase search_database;
		search_database.initialize(test_residues);
		TR.Trace << json(test_residues).dump(2) << std::endl;

		ndarray::Array<SearchReal, 3, 3> test_query_coordinates(5, 4, 3);
		test_query_coordinates.deep() = orient_array(test_residues)[ndarray::view(3,8)()()];

		SingleQueryExecutor query_executor(StructureSingleQuery(test_query_coordinates, .01));

		query_executor.execute(search_database);
		TR << "query summary: " << json(query_executor.query_stats) << std::endl;

		TS_ASSERT_EQUALS(query_executor.query_results.size(), 1);
		StructureSingleQueryResult result = query_executor.query_results.at(0);
		TS_ASSERT_EQUALS(result.fragment_start, 3);
		TS_ASSERT_DELTA(result.result_rmsd, 0, 1e-3);
	}

	void test_pair_query()
	{
		core::pose::PoseOP test_pose = core::import_pose::pose_from_file("protocols/indexed_structure_store/test_structure.pdb", core::import_pose::PDB_file);

		ndarray::Array<ResidueEntry, 1> test_residues = extract_residue_entries( *test_pose );

		StructureDatabase search_database;
		search_database.initialize(test_residues);
		TR.Trace << json(a_to_v(test_residues)).dump(2) << std::endl;

		ndarray::Array<SearchReal, 3, 3> test_query_coordinates_a(2, 4, 3);
		ndarray::Array<SearchReal, 3, 3> test_query_coordinates_b(2, 4, 3);

		test_query_coordinates_a.deep() = orient_array(test_residues)[ndarray::view(3,5)()()];
		test_query_coordinates_b.deep() = orient_array(test_residues)[ndarray::view(8,10)()()];

		PairQueryExecutor query_executor(StructurePairQuery(test_query_coordinates_a, test_query_coordinates_b, .1));

		query_executor.execute(search_database);
		TR << "query summary: " << json(query_executor.query_stats) << std::endl;

		TS_ASSERT_EQUALS(query_executor.query_results.size(), 1);
		StructurePairQueryResult result = query_executor.query_results.at(0);
		TS_ASSERT_EQUALS(result.fragment_a_start, 3);
		TS_ASSERT_EQUALS(result.fragment_b_start, 8);
		TS_ASSERT_DELTA(result.result_rmsd, 0, 1e-3);
	}


};
}
