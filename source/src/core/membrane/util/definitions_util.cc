// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       core/membrane/util/definitions_util.cc
///
/// @brief      Utility Functions for membrane definitions
///	@details    The membrane definition utility functions are designed to be invariants that
///             guarantee safe initialization, resetting and copying based on their
///             intended definition and should be used
///
/// @note       When adding a definition to definition.hh, you should
///			 	ALWAYS write an initialize funciton in this file based
///			  	on your definition specifications
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_util_types_util_cc
#define INCLUDED_core_membrane_util_types_util_cc

// Unit Headers
#include <core/membrane/util/definitions_util.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

namespace core {
namespace membrane {
namespace util {

/// @brief EmbedConfigInfo initialization
/// @param [none]
EmbedConfigInfoOP init_embedConfigInfo()
{
	// Create a new embed info object
	EmbedConfigInfoOP embedding = new EmbedConfigInfo;

	// Initialize Members of Embed Info
	embedding->center.assign(0, 0, 0);
	embedding->normal.assign(0, 0, 0);

	embedding->depth = 0;

	// Done!
	return embedding;
}

/// @brief EmbedSearchParams initialization
/// @param [none]
EmbedSearchParamsOP init_EmbedSearchParams()
{
	// Create new search info object
	EmbedSearchParamsOP embed_search = new EmbedSearchParams;
    
    // Initialize fields
	embed_search->center_search = false;
	embed_search->center_max_delta = 5;

	embed_search->normal_search = false;
	embed_search->normal_start_angle = 10;
	embed_search->normal_max_angle = 10;
	embed_search->normal_delta_angle = 10;
    
    embed_search->penalties = false;
    embed_search->no_interpolate_Mpair = false;

	// Return new object
	return embed_search;
}

} // util
} // membrane
} // core


#endif // INCLUDED_core_membrane_util_types_util_cc

