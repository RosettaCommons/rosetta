// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefFallbackConfiguration.hh
///
/// @brief      Fallback configuration for embed def - membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_hh
#define INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_hh

// Unit Headers
#include <core/membrane/io/EmbedDefFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Definition Fallback Configuration
/// @details Fallback configuration for embedding definitions if no resource specified in .xml
class EmbedDefFallbackConfiguration : public basic::resource_manager::FallbackConfiguration
{

public: // typedefs
	typedef basic::resource_manager::ResourceDescription ResourceDescription;

public: // functions

    /// @brief Constructor
	EmbedDefFallbackConfiguration();
    
    /// @brief Destructor
	virtual ~EmbedDefFallbackConfiguration();

    /// @brief Fallback to commandline option -in:file:embedfile
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

    /// @brief Return corresponding resource loader EmbedDefLoader
	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

    /// @brief Return locator id of commandline specified option
	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

    /// @brief Return resource options for embedding definitions
	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

    /// @brief Throw error message if required resource not specified (this is a required resource)
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:

    /// @brief Grab locator id from -in:file:embedfile option
	basic::resource_manager::LocatorID get_embedfile_from_options() const;

}; // class Embed Def Fallback Configuration

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_hh

