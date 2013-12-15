// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/func/GaussianFunc.cxxtest.hh
/// @brief  test suite for GaussianFunc function
/// @author James Thompson
/// @author modified Apr 23 2008 by Sergey Lyskov: rewriting to use UTracer

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/scoring/func/GaussianFunc.hh>
#include <core/scoring/func/GaussianFunc.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/types.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

//Auto Headers
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.GaussianFunc.cxxtest");

using namespace core;

class GaussianFuncTests : public CxxTest::TestSuite
{

public:
	GaussianFuncTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {

	}

	///////////////////////////////////////////////////////////////////////////////
	// ------------------------------------------ //
	/// @brief simple test minimization
	void test_gaussian_func()
	{
		using namespace core::scoring::func;

		// UTracer log file
		test::UTracer UT("core/scoring/func/GaussianFunc.u");

		GaussianFuncOP func( new GaussianFunc( 5.16, 1.5 ) );

		float const TOLERATED_ERROR = 0.001;
		const core::Real start = 2;
		const core::Real end   = 20;
		const core::Real res   = 0.5;

		UT.abs_tolerance( TOLERATED_ERROR ) << "\n";

		core::Size nsteps = core::Size( ( end - start ) / res );
		for ( core::Size i = 0; i < nsteps; ++i ) {
			core::Real r = start + (i * res);
			UT << "r=" << r << " func=" << func->func(r) << " dfunc=" << func->dfunc(r) << std::endl;
		}
	} // test_gaussian_func
};
