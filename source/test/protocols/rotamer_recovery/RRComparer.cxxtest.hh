// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RRComparer.cxxtest.hh
/// @brief  Comparison functions for rotamer_recovery
/// @author Matthew O'Meara (mattjomeara@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <util/pose_funcs.hh>

// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh>

// Project Headers
#include <test/core/init_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <ostream>
#include <sstream>
#include <vector>
#include <basic/Tracer.fwd.hh>



static basic::Tracer TR("protocols.rotamer_recovery.RRComparer.cxxtest");

class RRComparerTests : public CxxTest::TestSuite {

public:

	void
	setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -patch_selectors VIRTUAL_BB -override_rsd_type_limit" ); // kinemages below made with dun02 library
	}

	void test_RRComparerRotBins_main() {
		do_test_RRComparerRotBins_identity();

	}

	void
	do_test_RRComparerRotBins_identity(){
		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using protocols::rotamer_recovery::RRComparerRotBins;


		RRComparerRotBins rrc;

		Pose pose ( fullatom_pose_from_string( pdb_string_1ten() ) );

		Real score;
		bool recovered;

		for( Size i=1; i <= pose.total_residue(); ++i ){

			rrc.measure_rotamer_recovery(
				pose, pose, pose.residue(i), pose.residue(i), score, recovered);

			TS_ASSERT( score == 0 );
			TS_ASSERT( recovered );
		}

		rrc.measure_rotamer_recovery(
			pose, pose, pose.residue(11), pose.residue(2), score, recovered);
		//TR << "ASP3 vs ASP11: score: " << score << " recovered: " << recovered << endl;
		TS_ASSERT( score == 2 );
		TS_ASSERT( recovered == false );

		rrc.measure_rotamer_recovery(
			pose, pose, pose.residue(8), pose.residue(27), score, recovered);
		//TR << "GLU8 vs GLU27: score: " << score << " recovered: " << recovered << endl;
		TS_ASSERT( score == 3 );
		TS_ASSERT( recovered == false );

		rrc.measure_rotamer_recovery(
			pose, pose, pose.residue(24), pose.residue(41), score, recovered);
		//TR << "PRO24 vs PRO41: score: " << score << " recovered: " << recovered << endl;
		TS_ASSERT( score == 0 );
		TS_ASSERT( recovered == true );

	}

	void test_RRComparerAutomorphicRMSD_main() {
		do_test_RRComparerAutomorphicRMSD_identity();

	}

	void do_test_RRComparerAutomorphicRMSD_identity(){
		using core::Real;
		using core::Size;
		using core::pose::Pose;
		using protocols::rotamer_recovery::RRComparerAutomorphicRMSD;


		RRComparerAutomorphicRMSD rrc;

		Pose pose ( fullatom_pose_from_string( pdb_string_1ten() ) );

		Real score;
		bool recovered;

		for( Size i=1; i <= pose.total_residue(); ++i ){

			rrc.measure_rotamer_recovery(
				pose, pose, pose.residue(i), pose.residue(i), score, recovered);

			TS_ASSERT_DELTA( score, 0, .0001 );
			TS_ASSERT( recovered );
		}

	}

};
