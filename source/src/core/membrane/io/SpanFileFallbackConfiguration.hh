// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileFallbackConfiguration.hh
///
/// @brief      Fallback configuration for span file loader/options - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_memrbane_io_SpanFileFallbackConfiguration_hh
#define INCLUDED_core_memrbane_io_SpanFileFallbackConfiguration_hh

// Unit Headers
#include <core/membrane/io/SpanFileFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Fallback Configuration for Span Files
/// @details Fallback to options system if span file resource not specified in .xml
class SpanFileFallbackConfiguration : public basic::resource_manager::FallbackConfiguration
{

public: // typedefs

	typedef basic::resource_manager::ResourceDescription ResourceDescription;

public: // functions

    /// @brief Constructor
	SpanFileFallbackConfiguration();

    /// @brief Specify Fallback given resource description
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

    /// @brief Return Span File Loader Type
	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

    /// @brief Return locator ID for span file
	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

    /// @brief Get SpanFile Options class
	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

    /// @brief Return error message given no resource definiiton (required)
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:

    /// @brief Get span file locator id from cmd option
	basic::resource_manager::LocatorID get_span_filename_from_options() const;

}; // class SpanFileFallbackConfiguraiton

} // io
} // membrane
} // core

#endif // INCLUDED_core_memrbane_io_SpanFileFallbackConfiguration_hh

