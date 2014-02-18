// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsLoader.fwd.hh
///
/// @brief      Embedding Search Parameters loader class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsLoader_fwd_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsLoader_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {
            
/// @brief Embedding Search Params Loader Class
/// @details Load Embedding search params data as a resource
class EmbedSearchParamsLoader;
typedef utility::pointer::owning_ptr< EmbedSearchParamsLoader > EmbedSearchParamsLoaderOP;
typedef utility::pointer::owning_ptr< EmbedSearchParamsLoader const > EmbedSearchParamsLoaderCOP;
    
} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsLoader_fwd_hh

