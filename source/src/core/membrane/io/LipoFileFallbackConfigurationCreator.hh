// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileFallbackConfigurationCreator.hh
///
/// @brief      Fallback Configuration - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileFallbackConfigurationCreator_hh
#define INCLUDED_core_membrane_io_LipoFileFallbackConfigurationCreator_hh

// unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace membrane {
namespace io {

/// @brief Fallback Configuraiton for Lipo Files
/// @details Fallback to commandline options if no resource specified
class LipoFileFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
{

public:

    /// @brief Return fallback class
	virtual
	basic::resource_manager::FallbackConfigurationOP
	create_fallback_configuration() const;

    /// @brief Return lipo file fallback type
	virtual
	std::string
	resource_description() const;

};

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileFallbackConfigurationCreator_hh

