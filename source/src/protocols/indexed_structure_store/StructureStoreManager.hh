// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStoreManager.hh
/// @brief Manages access to structure store via name-specific backends.
/// @author Alex Ford (fordas@uw.edu)
//

#ifndef INCLUDED_protocols_indexed_structure_store_StructureStoreManager_hh
#define INCLUDED_protocols_indexed_structure_store_StructureStoreManager_hh

#include <mutex>
#include <deque>

// Project headers
#include <core/types.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.fwd.hh>
#include <protocols/indexed_structure_store/StructureStoreProvider.fwd.hh>
#include <protocols/indexed_structure_store/StructureStore.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <numeric/types.hh>
#include <iosfwd>
#include <map>
#include <tuple>
#include <utility/vector1.fwd.hh>

namespace protocols
{
namespace indexed_structure_store
{

// @brief Core database handle.
// Encapsulates reading Structure/Residue data from data store and manages retrieval on indices on store.
class StructureStoreManager : public utility::SingletonBase< StructureStoreManager >
{
public:
	friend class utility::SingletonBase< StructureStoreManager >;

	StructureStoreOP load_structure_store(std::string store_path);
	void write_structure_store(std::string store_path, StructureStore & store);

	void register_store_provider(core::SSize priority, std::string name, StructureStoreProviderOP backend);

private:
	StructureStoreManager();

	std::map<std::tuple<core::SSize, std::string>, StructureStoreProviderOP> store_providers;
	std::map<std::string, StructureStoreOP> store_cache;


	std::mutex cache_mutex;
};

}
}
#endif
