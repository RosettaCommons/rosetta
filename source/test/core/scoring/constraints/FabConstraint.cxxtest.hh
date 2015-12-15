// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/FabConstraint.hh.cxxtest.hh
/// @brief  test suite for FabConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <core/scoring/constraints/FabConstraint.hh>

#ifdef	SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <test/util/pose_funcs.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class FabConstraintTests : public CxxTest::TestSuite
{
public:

	void setUp() {
		core_init();
	}

	void test_serialize_FabConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::pose;
		using namespace core::id;

		Pose trpcage = create_trpcage_ideal_pose();
		FabConstraintOP fab_cst( new FabConstraint() ); // serialize this through a pointer to the base class
		utility::vector1< core::Size > r1, r2;
		r1.push_back( 3 ); r1.push_back( 4 );
		r2.push_back( 5 ); r2.push_back( 6 );
		fab_cst->setup_csts( trpcage, r1, r2, "AA" );
		// "AA" says that both the antibody and the antigen are the same chain;
		// this is an abuse of the constraint for testing purposes;

		// serialize a pointer to the base class
		ConstraintOP instance( fab_cst );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a FabConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< FabConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

