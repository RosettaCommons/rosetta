// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/H5StructureStoreBackend
/// @brief HDF5-backed protein structure store.
/// @author Alex Ford <fordas@uw.edu>
//

#ifndef INCLUDED_protocols_indexed_structure_store_H5StructureStoreBackend_hh
#define INCLUDED_protocols_indexed_structure_store_H5StructureStoreBackend_hh

#ifdef USEHDF5

// Utility Headers
#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/indexed_structure_store/H5StructureStoreBackend.fwd.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreProvider.hh>

#include <vector>
#include <string>

#include "H5Cpp.h"

namespace protocols
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from data store and manages retrieval on indices on store.
class H5StructureStoreBackend : public StructureStoreProvider
{
public:
	// Structure database contains structure data and structure geometry indices for a collection
	// of target structures. The data is stored in the following tables:
	//    '/structures' =
	//        [int structure_id, string structure_name]
	//    '/residues'   =
	//        [int structure_id, int residue_id, <residue identity data>, <residue torsion angles>]
	//
	H5StructureStoreBackend() = default;

  bool can_load(std::string store_path) override;
	bool can_write(std::string store_path) override;

	// @brief Retrieves structure and residue information from the backing store.
	StructureStoreOP load_store(std::string store_path) override;
	void write_store(std::string store_path, StructureStore & store) override;

	//H5 DataTypes for store datatypes.
	static H5::DataType StructureEntry_datatype();

	static H5::DataType ResidueEntry_datatype();
	static H5::DataType ResidueBackboneEntry_datatype();
	static H5::DataType ResidueSidechainEntry_datatype();
	static H5::DataType ResidueOrientEntry_datatype();
};

}
}

#endif
#endif
