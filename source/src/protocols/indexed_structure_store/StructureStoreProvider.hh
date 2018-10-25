// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStoreProvider.hh
/// @brief Abstract structure store backing interface.
/// @author Alex Ford (fordas@uw.edu)
//

#pragma once

// Utility Headers
#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/indexed_structure_store/StructureStoreProvider.fwd.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>

#include <vector>
#include <string>

namespace protocols
{
namespace indexed_structure_store
{

class StructureStoreProvider : public utility::pointer::ReferenceCount
{
public:
	virtual bool can_load(std::string store_path) = 0;
	virtual bool can_write(std::string store_path) = 0;

	// @brief Retrieves structure and residue information from the backing store.
	virtual StructureStoreOP load_store(std::string store_path) = 0;

	// @brief Write structure information to the backing store.
	virtual void write_store(std::string store_path, StructureStore & store) = 0;
};

}
}
