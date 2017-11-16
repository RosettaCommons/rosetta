// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/constraints/AntibodyConstraintTests.cxxtest.hh
/// @brief  Utility functions for testing residue selectors
/// @author Jared Adolf-Bryfogle

#ifndef INCLUDED_core_select_residue_selector_utilities_for_testing_HH
#define INCLUDED_core_select_residue_selector_utilities_for_testing_HH

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>


// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>
#include <iostream>



static basic::Tracer TR_util("core.selector.residue_selector.utilities_for_testing");




inline
void
compare_bool_vector(utility::vector1< bool > const & vec1, utility::vector1< bool > const & vec2){
	TS_ASSERT(vec1.size() == vec2.size());
	for ( core::Size i = 1; i <= vec1.size(); ++i ) {
		//

		//Debugging:
		if ( vec1[ i ] != vec2[ i ] ) {
			TR_util <<"NE "<< i << " " << utility::to_string(vec1[ i ]) << " " << utility::to_string(vec2[ i ]) << std::endl;
			for ( core::Size x = 1; x <= vec1.size(); ++x ) {
				TR_util << utility::to_string(vec1[ i ] != vec2[ i ]) << " " << x << " " << utility::to_string(vec1[ x ]) << " " << utility::to_string(vec2[ x ]) << std::endl;
			}
			TR_util.flush();

		}
		TS_ASSERT(vec1[ i ] == vec2[ i ]);

	}
}

#endif




