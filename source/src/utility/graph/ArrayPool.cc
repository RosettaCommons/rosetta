// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/ArrayPool.cc
/// @brief  ArrayPool class declaration and implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <utility/graph/ArrayPool.hh>

// Package Headers
#include <platform/types.hh>

// Utility headers
#include <utility/string_util.hh>

namespace utility {
namespace graph {

std::string
neg_space_element_allocation_error_message( platform::Size block_size, platform::Size neg_space_element_size )
{
	return "ERROR: new failed in ArrayPool when requesting neg space block of size " +
		utility::to_string( block_size ) + " (" +
		utility::to_string( block_size * neg_space_element_size ) +
		" bytes)";
}

std::string
block_allocation_error_message( platform::Size block_size, platform::Size array_size, platform::Size t_size )
{
	return "ERROR: new failed in ArrayPool when requesting Array block of size " +
		utility::to_string( block_size ) + "x" + utility::to_string( array_size ) + " (" +
		utility::to_string( block_size * array_size * t_size ) +
		" bytes)";
}

} //end namespace graph
} //end namespace utility
