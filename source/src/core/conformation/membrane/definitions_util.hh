// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       core/conformation/membrane/definitions_util.hh
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

#ifndef INCLUDED_core_conformation_membrane_definitions_util_hh
#define INCLUDED_core_conformation_membrane_definitions_util_hh

// Project Headers
#include <core/conformation/membrane/definitions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

namespace core {
namespace conformation {
namespace membrane {

//////// SpanningTopology Functions //////////

/// @brief EmbedConfigInfo initialization
/// @param [none]
EmbedConfigInfoOP init_embedConfigInfo();

/// @brief EmbedSearchParams initialization
/// @param [none]
EmbedSearchParamsOP init_EmbedSearchParams();

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_util_hh

