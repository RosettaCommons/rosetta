// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/BigBinConstraint.cxxtest.hh
/// @brief  test suite for BigBinConstraint function
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/constraints/BigBinConstraint.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/types.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/vector1.hh>
#include <iosfwd>
#include <string>

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class BigBinConstraintTests : public CxxTest::TestSuite {

public:
	BigBinConstraintTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_read_def() {
		using core::Size;
		using core::Real;
		using core::pose::Pose;
		using utility::vector1;
		using namespace core::scoring::constraints;

		vector1< Real > means, sdevs, weights;
		std::string const def( "BigBin 5 G 0.5" );
		std::istringstream input( def );

		BigBinConstraint bb_cst;

		Pose mypose;
		core::pose::make_pose_from_sequence(
			mypose,
			"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTL",
			"fa_standard"
		);

		core::scoring::func::FuncFactory ff; // dummy

		bb_cst.read_def( input, mypose, ff );

		TS_ASSERT( bb_cst.res() == 5 );
		TS_ASSERT( bb_cst.bin() == 'G' );
		TS_ASSERT_DELTA( bb_cst.sdev(), 0.5, 1e-5 );
	} // test_func

	void test_serialize_BigBinConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::id;

		// serialize this through a pointer to the base class
		ConstraintOP instance( new BigBinConstraint( AtomID(3,3), AtomID(1,4), AtomID(2,4), AtomID(3,4), AtomID(1,5), AtomID(2,5), 'E' ) );

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

		// make sure the deserialized base class pointer points to a BigBinConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< BigBinConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};
