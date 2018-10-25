// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStore.cc
/// @brief A POD-based database of residue-level protein structures.
/// @author Alex Ford (fordas@uw.edu)
//

#include <protocols/indexed_structure_store/StructureStore.hh>

#include <algorithm>
#include <vector>
#include <string>
#include <stdexcept>

namespace protocols
{
namespace indexed_structure_store
{

StructureStore::StructureStore()
{
}

void StructureStore::check_entries()
{
	for ( std::size_t i = 0; i < structure_entries.size() - 1; i++ ) {
		if ( structure_entries[i].id > structure_entries[i + 1].id ) {
			throw std::invalid_argument("structure_entries not read in sorted order");
		}

		if ( structure_entries[i].id == structure_entries[i + 1].id ) {
			throw std::invalid_argument("structure_entries not unique");
		}
	}

	for ( std::size_t i = 0; i < residue_entries.size() - 1; i++ ) {
		if ( ResidueEntry::comp_id(residue_entries[i + 1], residue_entries[i]) ) {
			throw std::invalid_argument("residue_entries not read in sorted order");
		}

		if ( !ResidueEntry::comp_id(residue_entries[i], residue_entries[i+1]) ) {
			throw std::invalid_argument("residue_entries not unique");
		}
	}
}

StructureEntry StructureStore::get_structure(uint32_t structure_id)
{
	return *std::lower_bound(
		structure_entries.begin(),
		structure_entries.end(),
		StructureEntry {structure_id, {0}},
		StructureEntry::comp_id
	);
}

ndarray::Array<const ResidueEntry, 1, 1> StructureStore::get_residues(uint32_t structure_id)
{
	ResidueEntry search_entry;
	search_entry.structure_id = structure_id;

	auto range = std::equal_range(
		residue_entries.begin(), residue_entries.end(),
		search_entry, [](const ResidueEntry a, const ResidueEntry b) { return a.structure_id < b.structure_id; }
	);

	int begin_i = std::get<0>(range) - residue_entries.begin();
	int end_i = std::get<1>(range) - residue_entries.begin();

	return residue_entries[ndarray::view(begin_i, end_i)];
}

void StructureStore::resize( uint32_t structures_size, uint32_t residues_size )
{
	ndarray::Array<StructureEntry, 1, 1> new_structure_entries = ndarray::allocate(structures_size);
	ndarray::Array<ResidueEntry, 1, 1> new_residue_entries = ndarray::allocate(residues_size);

	new_structure_entries[ndarray::view(0, structure_entries.getSize<0>())] = structure_entries;
	new_residue_entries[ndarray::view(0, residue_entries.getSize<0>())] = residue_entries;

	structure_entries = new_structure_entries;
	residue_entries = new_residue_entries;
}

}
}
