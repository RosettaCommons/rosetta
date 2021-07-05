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

// Package headers
#include <protocols/simple_filters/InterfaceHbondsFilter.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>


// Utility headers

// Basic headers
#include <basic/Tracer.hh>

// C++ headers

#include <core/init_util.hh> // AUTO IWYU For core_init

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
