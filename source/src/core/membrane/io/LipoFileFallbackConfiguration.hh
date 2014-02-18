// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileFallbackConfiguration.hh
///
/// @brief      Fallback Configuration - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_hh
#define INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_hh

// Unit headers
#include <core/membrane/io/LipoFileFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace core {
namespace membrane {
namespace io {
    
/// @brief Lipo File Fallback Configuration
/// @details Options cmd fallback if resource not specified
class LipoFileFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {

public:
	typedef basic::resource_manager::ResourceDescription ResourceDescription;

public:

    /// @brief Constructor
	LipoFileFallbackConfiguration();

    /// @brief Specify Option to Fallbck (-in:file:lipofile)
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

    /// @brief Get apporpriate resource loader (LipoFileLoader)
	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

    /// @brief Provide locator id for .lips4 file
	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

    /// @brief Get options class for lipo files
	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

    /// @brief throw error message if nothing is specified
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:

    /// @brief Grab lipo file id from opts
	basic::resource_manager::LocatorID get_lips_filename_from_options() const;

};

} // io
} // membrane
} // core


#endif // INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_hh

