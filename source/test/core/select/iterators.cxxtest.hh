// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/select/iterators.cxxtest.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>

#include <utility/iterators.hh>
#include <core/select/iterators.hh>

#include <test/core/init_util.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>


#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/kinematics/FoldTree.hh>

#include <basic/Tracer.hh>

#include <test/util/pose_funcs.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

static basic::Tracer TR("core.select.iterators.cxxtest");

using namespace utility;
using namespace core::select;

class PoseIteratorTests : public CxxTest::TestSuite {
public:
	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_enumerate1(){

		utility::vector1< bool > this_residue_has_been_sampled( 10, false );
		//Start out with each of these as false, then activate each one in the for loop

		using namespace core::pose;
		for ( core::Size const resid : enumerate1( 10 ) ) {
			this_residue_has_been_sampled[ resid ] = true;
		}

		for ( bool const marker : this_residue_has_been_sampled ) {
			TS_ASSERT( marker );
		}
	}

	void test_iterators(){
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "STEVENLEWIS", "fa_standard" );
		TS_ASSERT_EQUALS( pose.size(), 11 );

		utility::vector1< bool > this_residue_has_been_sampled( pose.size(), false );
		//Start out with each of these as false, then activate each one in the for loop
		using namespace core::pose;
		for ( core::Size const resid : resids( pose ) ) {
			this_residue_has_been_sampled[ resid ] = true;
		}

		for ( bool const marker : this_residue_has_been_sampled ) {
			TS_ASSERT( marker );
		}
	}

	void test_iterators_with_selectors(){
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "STEVENLEWIS", "fa_standard" );
		TS_ASSERT_EQUALS( pose.size(), 11 );

		utility::vector1< bool > this_residue_has_been_sampled( pose.size(), false );
		//Start out with each of these as false, then activate each one in the for loop

		using namespace core::select::residue_selector;
		ResidueIndexSelectorOP s1 =
			utility::pointer::make_shared< ResidueIndexSelector >( "1,2,3,4,5" );
		ResidueIndexSelectorOP s2 =
			utility::pointer::make_shared< ResidueIndexSelector >( "1,5,6,7" );
		//Only 1 and 5 should be selected

		using namespace core::pose;
		for ( core::Size const resid : resids( pose, {s1,s2} ) ) {
			this_residue_has_been_sampled[ resid ] = true;
		}

		for ( core::Size const index : indices1( this_residue_has_been_sampled ) ) {
			if ( index == 1 || index == 5 ) {
				TS_ASSERT( this_residue_has_been_sampled[ index ] );
			} else {
				TS_ASSERT( ! this_residue_has_been_sampled[ index ] );
			}
		}
	}

};

