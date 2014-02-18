// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileOptions.fwd.hh
///
/// @brief      Options for generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileOptions_fwd_hh
#define INCLUDED_core_membrane_io_SpanFileOptions_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {
    
    /// @brief Options for generating per-chain membrane spanning data
    /// @details Generate spanning topology object options set
    class SpanFileOptions;
    typedef utility::pointer::owning_ptr< SpanFileOptions > SpanFileOptionsOP;
    typedef utility::pointer::owning_ptr< SpanFileOptions const > SpanFileOptionsCOP;

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileOptions_fwd_hh
