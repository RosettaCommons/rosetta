// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/util/definitions.fwd.hh
///
/// @brief 		Definitions for membrane protein modeling data
/// @details 	Useful customd ata structures for implementing membrane protein
///				modeling framework
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_util_definitions_fwd_hh
#define INCLUDED_core_membrane_util_definitions_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace membrane {
namespace util {

/// @brief      Struct EmbedConfigInfo
/// @details    Stores data required for properly initializing membrane protein embedding
///             data either from user definition or external methods such as searc/score
struct EmbedConfigInfo;

typedef utility::pointer::owning_ptr< EmbedConfigInfo > EmbedConfigInfoOP;
typedef utility::pointer::owning_ptr< EmbedConfigInfo const > EmbedConfigInfoCOP;

/// @brief      Struct Embedding Search Info
/// @details	Stores required parameters for performing a search and score for embedding (protein)
///             in a membrane
struct EmbedSearchParams;

typedef utility::pointer::owning_ptr< EmbedSearchParams > EmbedSearchParamsOP;
typedef utility::pointer::owning_ptr< EmbedSearchParams const > EmbedSearchParamsCOP;

} // util
} // membrane
} // core

#endif // INCLUDED_core_membrane_util_definitions_fwd_hh

