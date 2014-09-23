// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#include <utility/exit.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <core/indexed_structure_store/StructureStoreManager.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>

// H5-based backend declarations are guarded by #ifdef USEHDF5 
#include <core/indexed_structure_store/H5FragmentStoreBackend.hh>
#include <core/indexed_structure_store/BinaryFragmentStoreBackend.hh>

namespace core
{
namespace indexed_structure_store
{

using namespace basic::options;

static thread_local basic::Tracer TR( "core.indexed_structure_store.StructureStoreManager" );

StructureStoreManager::StructureStoreManager() 
{
} 

/// @brief set initial value as no instance
StructureStoreManager* StructureStoreManager::instance_( 0 );

/// @brief static function to get the instance of ( pointer to) this singleton class
StructureStoreManager * StructureStoreManager::get_instance()
{
	if ( instance_ == 0 )
	{
		 instance_ = new StructureStoreManager();
	}
	return instance_;
}

FragmentLookupOP StructureStoreManager::load_fragment_lookup(std::string lookup_name)
{
	if( option[OptionKeys::indexed_structure_store::fragment_store])
	{
		std::string resolved_name = resolve_store_path(option[OptionKeys::indexed_structure_store::fragment_store]());
		std::string store_name = option[OptionKeys::indexed_structure_store::fragment_store]();

		if(resolved_name.empty())
		{
			utility_exit_with_message(
					"Unable to resolve store specified in indexed_structure_store::fragment_store: " + store_name );
		}

		return load_fragment_lookup(lookup_name, resolved_name);
	}
	else
	{
			utility_exit_with_message("Specify indexed_structure_store::fragment_store target file.");
	}
}

FragmentLookupOP StructureStoreManager::load_fragment_lookup(std::string lookup_name, std::string store_path)
{
	FragmentStoreOP target_store(NULL);

	std::string resolved_path = resolve_store_path(store_path);

	if(resolved_path.empty())
	{
			utility_exit_with_message("Unable to resolve specified store: " + store_path);
	}

	if(utility::file::file_extension(resolved_path) == "h5")
	{
#ifdef USEHDF5
		H5FragmentStoreBackend backend(store_path);
		target_store = backend.get_fragment_store(lookup_name);
#else
		utility_exit_with_message("StructureStoreManager::load_fragment_lookup without HDF5 support, unable to load lookup: " + lookup_name + " resolved from: " + store_path);
#endif
	}
	else
	{
		BinaryFragmentStoreBackend backend(resolved_path);
		target_store = backend.get_fragment_store(lookup_name);
	}

	FragmentLookupOP lookup( new FragmentLookup(target_store) );
	return lookup;
}

std::string StructureStoreManager::resolve_store_path(std::string target_path)
{
	std::string basename = utility::file::file_basename(target_path);

#ifdef USEHDF5
	// Resolve h5 stores first, given either the basename or full path
	if(utility::file::file_exists(basename + ".h5"))
	{
		return basename + ".h5";
	}
#endif

	// Fall back to target or basename
	if(utility::file::file_exists(target_path))
	{
		return target_path;
	}

	if(utility::file::file_exists(basename))
	{
		return basename;
	}

	return "";
}

}
}
