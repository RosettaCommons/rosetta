// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  Test suite for SEWING Assembly class
/// @author Tim Jacobs (TimJacobs2@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <devel/sewing/DisembodiedAssembly.hh>
#include <devel/sewing/Hasher.hh>
#include <devel/sewing/SewGraph.hh>
#include <devel/sewing/util.hh>
#include <devel/sewing/io.hh>
#include <test/core/init_util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility headers
#include <utility/sql_database/DatabaseSessionManager.hh>

/// Project headers

//Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel headers
#include <devel/init.hh>

// --------------- Test Class --------------- //
using namespace devel::sewing;

class AssemblyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();

        //Generate a graph to be used by all the tests
        utility::vector1<BasisPair> alignment_pairs = read_hashing_scores_from_file("devel/sewing/inputs/test.scores");
        //std::cout << "Number of alignment pairs: " << alignment_pairs.size() << std::endl;
        std::map<core::Size, Model> models = read_model_file("devel/sewing/inputs/test.models");
        continuous_graph_ = new SewGraph(models, 1);
        //continuous_graph_->add_all_model_edges_from_binary(edges)
	}

	// Shared finalization goes here.
	void tearDown() {
	}

    void test_graph_generation() {
        TS_ASSERT_EQUALS(continuous_graph_->num_nodes(), 2);
        TS_ASSERT_EQUALS(continuous_graph_->num_edges(), 1);
        TS_ASSERT_EQUALS(discontinuous_graph_->num_nodes(), 2);
        TS_ASSERT_EQUALS(discontinuous_graph_->num_edges(), 1);
    }

    void test_continuous_append() {
    }

    void test_discontinuous_append() {
    }

private:
    SewGraphOP continuous_graph_;
    SewGraphOP discontinuous_graph_;
};
