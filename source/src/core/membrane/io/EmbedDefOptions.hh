// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefOptions.hh
///
/// @brief      Options for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefOptions_hh
#define INCLUDED_core_membrane_io_EmbedDefOptions_hh

// Unit Headers
#include <core/membrane/io/EmbedDefOptions.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceOptions.hh>

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <cstdlib>
#include <string>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Definition Options
/// @details Options class for embedding definitions for protein chains
class EmbedDefOptions : public basic::resource_manager::ResourceOptions
{

public:

    /// @brief Constructor
	EmbedDefOptions();
    
    /// @brief Destructor
	virtual ~EmbedDefOptions();

    /// @brief Parse options from .xml file
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
		);

    /// @brief Return options class type
	virtual
	std::string
	type() const;


}; // class Embed Definition options

} // io
} // membrane
} // core


#endif // INCLUDED_core_membrane_io_EmbedDefOptions_hh


