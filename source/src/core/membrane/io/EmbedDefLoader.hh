// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefLoader.hh
///
/// @brief      Loader class for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefLoader_hh
#define INCLUDED_core_membrane_io_EmbedDefLoader_hh

// Unit Header
#include <core/membrane/io/EmbedDefLoader.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <istream>

namespace core {
namespace membrane {
namespace io {

/// @brief Loader for Embedding Definitions
/// @details Resource Loader for membrane protein chain embeddings
class EmbedDefLoader : public basic::resource_manager::ResourceLoader {

public:

	/// @brief Constructor
	EmbedDefLoader();

	/// @brief Destructor
	virtual ~EmbedDefLoader();

	/// @brief Returns an OP to embedding search options
	virtual
	utility::pointer::ReferenceCountOP
	create_resource(
		basic::resource_manager::ResourceOptions const & options,
		basic::resource_manager::LocatorID const & locator_id,
		std::istream & istream
		) const;

	/// @brief Default options
	virtual
	basic::resource_manager::ResourceOptionsOP
	default_options() const;

}; // class EmbedDefLoader

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefLoader_hh

