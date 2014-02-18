// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsFallbackConfiguration.fwd.hh
///
/// @brief      Embedding Search Parameters Fallback Configuration class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguraiton_fwd_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguraiton_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {
    
/// @brief Fallback Configuration for Embedding Search Parameters
/// @details Fallback if no resource for embedding search parameters specified
class EmbedSearchParamsFallbackConfiguraiton;
typedef utility::pointer::owning_ptr< EmbedSearchParamsFallbackConfiguraiton > EmbedSearchParamsFallbackConfiguraitonOP;
typedef utility::pointer::owning_ptr< EmbedSearchParamsFallbackConfiguraiton const > EmbedSearchParamsFallbackConfiguraitonCOP;
            
} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguraiton_fwd_hh


