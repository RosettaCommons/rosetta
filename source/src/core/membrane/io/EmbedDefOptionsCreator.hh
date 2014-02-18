// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefOptionsCreator.hh
///
/// @brief      Options Creator class for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefOptionsCreator_hh
#define INCLUDED_core_membrane_io_EmbedDefOptionsCreator_hh

// Package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Options Creator for Embedding Definitions
/// @details Options creator for embedding definitions - registered with core init
class EmbedDefOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{

public:

    /// @brief Return options type - embed def
	virtual
	std::string
	options_type() const;

    /// @brief Return new options class - embed def
	virtual
	basic::resource_manager::ResourceOptionsOP
	create_options() const;

}; // class EmbedDef options creator


} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefOptionsCreator_hh

