// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/indexed_structure_store/FragmentLookup.cxxtest.hh
/// @brief Test suite for core/indexed_structure_store/FragmentLookup
/// @author Alex Ford (fordas@uw.edu)
//
// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>
#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <core/indexed_structure_store/StructureStoreManager.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>

#include <core/indexed_structure_store/BinaryFragmentStoreBackend.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <vector>
#include <iterator>


static basic::Tracer TR("core.indexed_structure_store.FragmentLookup.cxxtest");

namespace
{
class FragmentLookupTests : public CxxTest::TestSuite {

public:
	void setUp()
	{
		core_init_with_additional_options( "" );
	}

	void test_read_store()
	{
		using namespace core::indexed_structure_store;
		core::pose::PoseOP test_pose = core::import_pose::pose_from_file("core/indexed_structure_store/test_structure.pdb", core::import_pose::PDB_file);

		FragmentLookupOP full_lookup = StructureStoreManager::get_instance()->load_fragment_lookup(
			"full_test_fragments", "core/indexed_structure_store/test_store");
		FragmentLookupOP partial_lookup = StructureStoreManager::get_instance()->load_fragment_lookup(
			"partial_test_fragments", "core/indexed_structure_store/test_store");

		std::vector<FragmentLookupResult> lookup_result;
		std::vector<numeric::Size> lookup_residue;

		lookup_result.resize(0);
		lookup_residue.resize(0);
		full_lookup->lookup_pose_fragments(
			*test_pose,
			std::back_inserter(lookup_result),
			std::back_inserter(lookup_residue));

		TS_ASSERT_EQUALS(
			lookup_result.size(),
			test_pose->size() - full_lookup->fragment_specification().fragment_length + 1);

		for ( numeric::Size i = 0; i < lookup_result.size(); i++ ) {
			TS_ASSERT_EQUALS(lookup_residue[i], i + 1);
			TS_ASSERT(lookup_result[i].found_match);
		}

		lookup_result.resize(0);
		lookup_residue.resize(0);
		partial_lookup->lookup_pose_fragments(
			*test_pose,
			std::back_inserter(lookup_result),
			std::back_inserter(lookup_residue));

		TS_ASSERT_EQUALS(
			lookup_result.size(),
			test_pose->size() - full_lookup->fragment_specification().fragment_length + 1);

		for ( numeric::Size i = 0; i < lookup_result.size() - 1; i++ ) {
			TS_ASSERT_EQUALS(lookup_residue[i], i + 1);
			TS_ASSERT(lookup_result[i].found_match);
		}

		numeric::Size f_i = lookup_result.size() - 1;
		TS_ASSERT_EQUALS(lookup_residue[f_i], f_i + 1);
		TS_ASSERT( !lookup_result[f_i].found_match );
	}

	void test_binary_store()
	{
		using namespace core::indexed_structure_store;
		core::pose::PoseOP test_pose = core::import_pose::pose_from_file("core/indexed_structure_store/test_structure.pdb", core::import_pose::PDB_file);

		FragmentSpecification test_spec;
		test_spec.fragment_atoms.push_back("N");
		test_spec.fragment_atoms.push_back("CA");
		test_spec.fragment_atoms.push_back("C");
		test_spec.fragment_length = 5;

		BinaryFragmentStoreBackend backend("core/indexed_structure_store/test_store");
		FragmentStoreOP full_lookup = backend.get_fragment_store("full_test_fragments");

		TS_ASSERT_EQUALS(full_lookup->fragment_specification.fragment_length, test_spec.fragment_length);
		TS_ASSERT_EQUALS(full_lookup->fragment_specification.fragment_atoms.size(), test_spec.fragment_atoms.size());
		TS_ASSERT_EQUALS(full_lookup->fragment_specification.fragment_atoms[0], test_spec.fragment_atoms[0]);
		TS_ASSERT_EQUALS(full_lookup->fragment_specification.fragment_atoms[1], test_spec.fragment_atoms[1]);
		TS_ASSERT_EQUALS(full_lookup->fragment_specification.fragment_atoms[2], test_spec.fragment_atoms[2]);

		numeric::Size full_expected_fragments = test_pose->size() - test_spec.fragment_length + 1;
		TS_ASSERT_EQUALS(full_lookup->fragment_threshold_distances.size(), full_expected_fragments);
		TS_ASSERT_EQUALS(full_lookup->fragment_coordinates.size(), full_expected_fragments * test_spec.coordinates_per_fragment());

		FragmentStoreOP partial_lookup = backend.get_fragment_store("partial_test_fragments");

		TS_ASSERT_EQUALS(partial_lookup->fragment_specification.fragment_length, test_spec.fragment_length);
		TS_ASSERT_EQUALS(partial_lookup->fragment_specification.fragment_atoms.size(), test_spec.fragment_atoms.size());
		TS_ASSERT_EQUALS(partial_lookup->fragment_specification.fragment_atoms[0], test_spec.fragment_atoms[0]);
		TS_ASSERT_EQUALS(partial_lookup->fragment_specification.fragment_atoms[1], test_spec.fragment_atoms[1]);
		TS_ASSERT_EQUALS(partial_lookup->fragment_specification.fragment_atoms[2], test_spec.fragment_atoms[2]);

		numeric::Size partial_expected_fragments = test_pose->size() - test_spec.fragment_length;
		TS_ASSERT_EQUALS(partial_lookup->fragment_threshold_distances.size(), partial_expected_fragments);
		TS_ASSERT_EQUALS(partial_lookup->fragment_coordinates.size(), partial_expected_fragments * test_spec.coordinates_per_fragment());
	}
};
}
