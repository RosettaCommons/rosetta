// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/PyReturnValuePolicyTest.hh
/// @brief A few functions test how PyRosetta handle boost ReturnValuePolicy
/// @author Sergey Lyskov


namespace protocols {
namespace toolbox {

// test for out-of-bounds access handling in Python
void out_of_bounds_memory_access()
{
	int *p = new int[1];
	for ( ; ; ++p ) *p=1;
	delete [] p;
}


} //toolbox
} //protocols


#include <protocols/toolbox/PyReturnValuePolicyTest.hh>
