// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceHbondsFilter.cxxtest.hh
/// @brief test suite for protocols::simple_filters::InterfaceHbondsFilter
/// @author Longxing Cao (longxing@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/simple_filters/InterfaceHbondsFilter.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rosetta_scripts/XmlObjects.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace protocols::simple_filters;

static basic::Tracer TR("protocols.simple_filters.InterfaceHbondsFilter");


class InterfaceHbondsFilterTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		core::import_pose::pose_from_file( pose_dimer_ab, "protocols/simple_filters/interface_hbonds_input.pdb" );
	}

	void test_InterfaceHbondsFilter() {
		TR << "test_InterfaceHbondsFilter_from_constructor_not_defaults" << std::endl;

		InterfaceHbondsFilterOP filter( new InterfaceHbondsFilter() );

		TS_ASSERT( 3.0 == filter->compute( pose_dimer_ab ) );

		filter->set_salt_bridge_mode( true );

		TS_ASSERT( 1.0 == filter->compute( pose_dimer_ab ) );

	}

private:
	core::pose::Pose pose_dimer_ab;

};
