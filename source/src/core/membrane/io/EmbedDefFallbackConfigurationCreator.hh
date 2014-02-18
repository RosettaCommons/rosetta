// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefFallbackConfigurationCreator.hh
///
/// @brief      Fallback configuration creator for embed def - membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefFallbackConfigurationCreator_hh
#define INCLUDED_core_membrane_io_EmbedDefFallbackConfigurationCreator_hh

// Unit Headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// Package Headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Definition Fallback Config Creator
/// @details Register embedding definition fallback with core init
class EmbedDefFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator {

public:
    
    /// @brief Specify fallback configuration class
	virtual
	basic::resource_manager::FallbackConfigurationOP
	create_fallback_configuration() const;

    /// @brief Return corresponting resource type - embed def
	virtual
	std::string resource_description() const;

};

} // io
} // membrane
} // core


#endif // INCLUDED_core_membrane_io_EmbedDefFallbackConfigurationCreator_hh

