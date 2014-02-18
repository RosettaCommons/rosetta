// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//// @file       EmbedDefLoaderreator.hh
///
/// @brief      Loader Creator class for loading membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefLoaderCreator_hh
#define INCLUDED_core_membrane_io_EmbedDefLoaderCreator_hh

// Unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Creator Class for the Embedding Definition Loader
/// @details Reigster the embedding definition loader class with core init
class EmbedDefLoaderCreator : public basic::resource_manager::ResourceLoaderCreator {
public:

    /// @brief Create embed def resource loader
	virtual
	basic::resource_manager::ResourceLoaderOP
	create_resource_loader() const;

    /// @brief Return laoder type embed def
	virtual
	std::string
	loader_type() const;

}; // class EmbedDefLoaderCreator

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefLoaderCreator_hh

