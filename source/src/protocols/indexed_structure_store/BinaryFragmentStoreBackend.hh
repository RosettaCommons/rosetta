// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#ifndef INCLUDED_protocols_indexed_structure_store_BinaryFragmentStoreBackend_hh
#define INCLUDED_protocols_indexed_structure_store_BinaryFragmentStoreBackend_hh


// Utility Headers
#include <vector>
#include <string>

#include <platform/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/indexed_structure_store/FragmentStoreProvider.hh>
#include <protocols/indexed_structure_store/BinaryFragmentStoreBackend.fwd.hh>
#include <protocols/indexed_structure_store/FragmentStore.fwd.hh>


namespace protocols
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from data store and manages retrieval on indices on store.
class BinaryFragmentStoreBackend : public FragmentStoreProvider
{
public:
	// Structure database contains structure data and structure geometry indices for a collection
	// of target structures. The data is stored in the following tables:
	//    TODO alexford document
	//
	// Opens file handle
	BinaryFragmentStoreBackend(std::string target_path);
	BinaryFragmentStoreBackend() = default;

	// @brief Retrieves fragment store from backend.
	FragmentStoreOP get_fragment_store(std::string store_name) override;

	void set_target_filename(std::string target_filename) override;
	void append_to_fragment_store(FragmentStoreOP fragment_store, std::string store_name, std::string group_field, std::string group_type) override;

private:
	std::string target_path_;
};

}
}

#endif
