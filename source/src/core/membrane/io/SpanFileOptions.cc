// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileOptions.cc
///
/// @brief      Options for generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileOptions_cc
#define INCLUDED_core_membrane_io_SpanFileOptions_cc

// Unit Headers
#include <core/membrane/io/SpanFileOptions.hh>
#include <core/membrane/io/SpanFileOptionsCreator.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace io {
    
    /// @brief Constructor
    SpanFileOptions::SpanFileOptions() : basic::resource_manager::ResourceOptions() {}
    
    /// @brief Destructor
    SpanFileOptions::~SpanFileOptions() {}

    /// @brief Parse Options from .xml Resource file
    void
    SpanFileOptions::parse_my_tag(
                                  utility::tag::TagCOP
                                  )
    {}

    /// @brief Return options class type
    std::string
    SpanFileOptions::type() const
    {
        return "SpanFileOptions";
    }

    /// @brief Return options class type to registrator
    std::string
    SpanFileOptionsCreator::options_type() const { return "SpanFileOptions"; }

    /// @brief Return options class to registrator
    basic::resource_manager::ResourceOptionsOP
    SpanFileOptionsCreator::create_options() const { return new SpanFileOptions; }

} // io
} // memrbane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileOptions_cc
