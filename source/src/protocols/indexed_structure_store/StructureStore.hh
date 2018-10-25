// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStore.hh
/// @brief A POD-based database of residue-level protein structures.
/// @author Alex Ford (fordas@uw.edu)
//
#ifndef INCLUDED_protocols_indexed_structure_store_StructureStore_hh
#define INCLUDED_protocols_indexed_structure_store_StructureStore_hh

// Utility Headers
#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/indexed_structure_store/StructureStore.fwd.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>

#include <ndarray.h>
#include <string>

namespace protocols
{
namespace indexed_structure_store
{

// @brief Core database handle.
class StructureStore : public utility::pointer::ReferenceCount
{
public:
	// @brief Basic structure store, holds a collection of structure and associated residue entries.
	StructureStore();

	// Retrieves structure entry for the given structure id.
	StructureEntry get_structure(uint32_t structure_id);

	// Retrieves array slice of residue entries.
	ndarray::Array<const ResidueEntry, 1, 1> get_residues(uint32_t structure_id);

	//@brief Validate entries by checking sort order.
	void check_entries();

	void resize(uint32_t structure_entries, uint32_t residue_entries);

	ndarray::Array<StructureEntry, 1, 1> structure_entries;
	ndarray::Array<ResidueEntry, 1, 1> residue_entries;
};

}
}
#endif
