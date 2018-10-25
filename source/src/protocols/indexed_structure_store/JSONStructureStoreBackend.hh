// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/JSONStructureStoreBackend
/// @brief JSON-backed protein structure store.
/// @author Alex Ford <fordas@uw.edu>
//

#pragma once

// Utility Headers
#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreProvider.hh>


#include <vector>
#include <string>

namespace protocols
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from concatenated JSON/Msgpack files.
class JSONStructureStoreBackend : public StructureStoreProvider
{
public:
	// @brief Retrieves structure and residue information from the backing store.
	virtual bool can_load(std::string store_path);
	virtual bool can_write(std::string store_path);
	// @brief Retrieves structure and residue information from the backing store.
	virtual StructureStoreOP load_store(std::string store_path);
	virtual void write_store(std::string store_path, StructureStore & store);
};

}
}
