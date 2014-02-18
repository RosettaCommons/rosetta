// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileFallbackConfigurationCreator.hh
///
/// @brief      Fallback configuration creator for span file loader/options - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileFallbackConfigurationCreator_hh
#define INCLUDED_core_membrane_io_SpanFileFallbackConfigurationCreator_hh

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace core {
namespace membrane {
namespace io {
    
/// @brief Register Span File Fallback
/// @details Provide class and description to core init
class SpanFileFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
{
public:
    
    /// @brief Return fallback class
	virtual
	basic::resource_manager::FallbackConfigurationOP
	create_fallback_configuration() const;

    /// @brief Return fallback class type
	virtual
	std::string resource_description() const;

};

} // io
} // membrane
} // core
#endif // INCLUDED_core_membrane_io_SpanFileFallbackConfigurationCreator_hh

