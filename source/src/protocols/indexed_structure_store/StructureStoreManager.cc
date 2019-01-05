// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStoreManager.cc
/// @brief Manages access to structure store via name-specific backends.
/// @author Alex Ford (fordas@uw.edu)
//

#include <utility/exit.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <protocols/indexed_structure_store/StructureStoreManager.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>

#include <protocols/indexed_structure_store/StructureStoreProvider.hh>
// H5-based backend declarations are guarded by #ifdef USEHDF5
#include <protocols/indexed_structure_store/H5StructureStoreBackend.hh>
#include <protocols/indexed_structure_store/JSONStructureStoreBackend.hh>
#include <protocols/indexed_structure_store/DirStructureStoreBackend.hh>
#include <protocols/indexed_structure_store/SilentStructureStoreBackend.hh>

#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace indexed_structure_store {

static basic::Tracer TR( "core.indexed_structure_store.StructureStoreManager" );

using namespace basic::options;
using utility::vector1;

StructureStoreManager::StructureStoreManager() : store_providers()
{
	typedef std::tuple<core::SSize, std::string> Key;
#ifdef USEHDF5
	store_providers[Key(0, "hdf5")] = StructureStoreProviderOP(new H5StructureStoreBackend());
#endif
	store_providers[Key( 5, "json" )] = utility::pointer::make_shared< JSONStructureStoreBackend >();
	store_providers[Key(10, "silent")] = utility::pointer::make_shared< SilentStructureStoreBackend >();
	store_providers[Key(10, "dir")] = utility::pointer::make_shared< DirStructureStoreBackend >();
}

StructureStoreOP StructureStoreManager::load_structure_store(std::string store_path)
{
	std::lock_guard<std::mutex> cache_lock(cache_mutex);

	TR << "Attempting to load store_path: " << store_path << std::endl;

	if ( store_cache.count(store_path) > 0 ) {
		try {

			TR << "Trying cached store." << std::endl;
			return store_cache[store_path];
		} catch (std::bad_weak_ptr &) {

			TR.Debug << "Cached store expired." << std::endl;
			store_cache.erase(store_path);
		}
	}

	for ( auto & provider : store_providers ) {
		auto & prio = std::get<0>(provider.first);
		auto & name = std::get<1>(provider.first);
		auto backend = provider.second;
		TR.Debug << "Checking backend: " << name << " prio: " << prio << std::endl;
		if ( backend->can_load(store_path) ) {
			TR << "Loading store: " << store_path  << " with backend: " << name << std::endl;
			StructureStoreOP store = backend->load_store(store_path);

			TR.Debug << "Caching store: " << store_path  << std::endl;
			store_cache[store_path] = store;

			return store;
		}
	}

	if ( utility::file::file_extension(store_path) == "h5" ) {
#ifndef USEHDF5
		utility_exit_with_message("StructureStoreManager::load_structure_store without HDF5 support, unable to load: " + store_path);
#endif
	}

	utility_exit_with_message("Unable to load specified store: " + store_path);
}

void StructureStoreManager::write_structure_store(std::string store_path, StructureStore & store) {
	for ( auto & provider : store_providers ) {
		auto & prio = std::get<0>(provider.first);
		auto & name = std::get<1>(provider.first);
		auto backend = provider.second;

		TR.Debug << "Checking backend: " << name << " prio: " << prio << std::endl;
		if ( backend->can_write(store_path) ) {
			TR << "Writing store: " << store_path  << " with backend: " << name << std::endl;
			backend->write_store(store_path, store);

			TR.Debug << "Clearing store: " << store_path  << std::endl;
			store_cache.erase(store_path);
			return;
		}
	}

	utility_exit_with_message("Unable to write specified store: " + store_path);
}

void StructureStoreManager::register_store_provider(
	core::SSize priority,
	std::string name,
	StructureStoreProviderOP backend)
{
	std::lock_guard<std::mutex> cache_lock(cache_mutex);

	for ( auto & provider : store_providers ) {
		auto & existing_name = std::get<1>(provider.first);
		if ( name == existing_name ) {
			utility_exit_with_message("Attempting to register duplicate provider name in StructureStoreManager.");
		}

		store_providers[std::make_tuple(priority, name)] = backend;
	}
}
}
}
