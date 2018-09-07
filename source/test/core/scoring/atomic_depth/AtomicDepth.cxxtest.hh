// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/atomic_depth/AtomicDepth.cxxtest.hh
/// @brief test suite for core::scoring::atomic_depth::AtomicDepth
/// @author Brian Coventry (bcov@uw.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/atomic_depth/util.hh>

// Project headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.tmpl.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

using namespace core;
using namespace core::scoring::atomic_depth;

static basic::Tracer TR("core.scoring.atomic_depth.AtomicDepthTests");


class AtomicDepthTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		core::import_pose::pose_from_file( pose_, "core/scoring/membrane/4O79_A.pdb" );
	}

	void check_depth( Real expected, Real actual ) {
		TR << "Expected: " << expected << " Calculated: " << actual << std::endl;

		Real delta = std::abs( expected - actual );
		TS_ASSERT( delta < 0.5f );
	}

	/// @brief Test probe radius 1.4
	void test_AtomicDepth_depth_test_14() {

		core::pose::Pose pose = pose_;

		AtomicDepth depth( pose_, 1.4, false, 0.25 );


		chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();

		utility::vector1<conformation::Atom> atoms = utility::vector1<conformation::Atom> {
			pose.residue(84).atom("CD1"),
			pose.residue(200).atom("CD2")
			};

		utility::vector1<Real> depths = depth.calcdepth( atoms, type_set );

		check_depth( 2.0, depths[1] );
		check_depth( 2.0, depths[2] );
	}

	/// @brief Test probe radius 4.0
	void test_AtomicDepth_depth_test_40() {

		core::pose::Pose pose = pose_;

		AtomicDepth depth( pose_, 4.0, false, 0.25 );


		chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();

		utility::vector1<conformation::Atom> atoms = utility::vector1<conformation::Atom> {
			pose.residue(84).atom("CD1"),
			pose.residue(200).atom("CD2")
			};

		utility::vector1<Real> depths = depth.calcdepth( atoms, type_set );

		check_depth( 2.0, depths[1] );
		check_depth( 8.0, depths[2] );
	}

	/// @brief Test probe radius 4.0
	void test_AtomicDepth_depth_test_polyleu_40() {

		core::pose::Pose pose = pose_;

		AtomicDepth depth( pose_, 4.0, true, 0.25 );


		chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();

		utility::vector1<conformation::Atom> atoms = utility::vector1<conformation::Atom> {
			pose.residue(84).atom("CD1"),
			pose.residue(200).atom("CD2")
			};

		utility::vector1<Real> depths = depth.calcdepth( atoms, type_set );

		check_depth( 2.0, depths[1] );
		check_depth( 7.75, depths[2] );
	}

	/// @brief Test probe radius 4.0
	void test_AtomicDepth_depth_test_oob() {

		core::pose::Pose pose = pose_;

		AtomicDepth depth( pose_, 4.0, false, 0.25 );


		chemical::AtomTypeSet const & type_set = pose.residue(1).type().atom_type_set();

		utility::vector1<conformation::Atom> atoms = utility::vector1<conformation::Atom> {
			pose.residue(200).atom("CD2"),
			};
		numeric::xyzVector<Real> oob( 10000, 10000, 10000);
		atoms[1].xyz( oob );

		utility::vector1<Real> depths = depth.calcdepth( atoms, type_set );

		check_depth( 2.0, depths[1] );
	}

	/// @brief Test probe radius 4.0
	void test_AtomicDepth_depth_test_util_40() {

		core::pose::Pose pose = pose_;
		AtomicDepthOP depth = nullptr;
		core::id::AtomID_Map< core::Real > depths;

		Size atno1 = pose.residue(84).atom_index("CD1");
		Size atno2 = pose.residue(200).atom_index("CD2");

		core::id::AtomID_Map< bool > depth_atoms;
		pose::initialize_atomid_map( depth_atoms, pose.conformation(), false );
		depth_atoms( 84, atno1 ) = true;
		depth_atoms( 200, atno2 ) = true;


		// Make sure all 4 functions work

		depths = atomic_depth( pose, 4.0 );

		check_depth( 2.0, depths(84, atno1) );
		check_depth( 8.0, depths(200, atno2) );

		depth = nullptr;
		depths = atomic_depth( pose, depth, 4.0 );

		check_depth( 2.0, depths(84, atno1) );
		check_depth( 8.0, depths(200, atno2) );
		TS_ASSERT( depth );

		depths = atomic_depth( pose, depth_atoms, 4.0 );

		check_depth( 2.0, depths(84, atno1) );
		check_depth( 8.0, depths(200, atno2) );

		depth = nullptr;
		depths = atomic_depth( pose, depth_atoms, depth, 4.0 );

		check_depth( 2.0, depths(84, atno1) );
		check_depth( 8.0, depths(200, atno2) );
		TS_ASSERT( depth );

		// One final check to make sure that depth caching works

		depths = atomic_depth( pose, depth_atoms, depth, 0.1 );

		check_depth( 2.0, depths(84, atno1) );
		check_depth( 8.0, depths(200, atno2) );

	}


private:
	core::pose::Pose pose_;

};
