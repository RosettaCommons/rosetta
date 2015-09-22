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

// Singleton instance and mutex static data members
namespace utility {

using core::indexed_structure_store::StructureStoreManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< StructureStoreManager >::singleton_mutex_{};
template <> std::atomic< StructureStoreManager * > utility::SingletonBase< StructureStoreManager >::instance_( 0 );
#else
template <> StructureStoreManager * utility::SingletonBase< StructureStoreManager >::instance_( 0 );
#endif

}

namespace core {
namespace indexed_structure_store {

static THREAD_LOCAL basic::Tracer TR( "core.indexed_structure_store.StructureStoreManager" );

using namespace basic::options;

StructureStoreManager::StructureStoreManager()
{
}

FragmentLookupOP StructureStoreManager::load_fragment_lookup(std::string lookup_name)
{
	if ( option[OptionKeys::indexed_structure_store::fragment_store] ) {
		std::string resolved_name = resolve_store_path(option[OptionKeys::indexed_structure_store::fragment_store]());
		std::string store_name = option[OptionKeys::indexed_structure_store::fragment_store]();

		if ( resolved_name.empty() ) {
			utility_exit_with_message(
				"Unable to resolve store specified in indexed_structure_store::fragment_store: " + store_name );
		}

		return load_fragment_lookup(lookup_name, resolved_name);
	} else {
		utility_exit_with_message("Specify indexed_structure_store::fragment_store target file.");
	}
}

FragmentLookupOP StructureStoreManager::load_fragment_lookup(std::string lookup_name, std::string store_path)
{
	FragmentStoreOP target_store(NULL);

	std::string resolved_path = resolve_store_path(store_path);

	if ( resolved_path.empty() ) {
		utility_exit_with_message("Unable to resolve specified store: " + store_path);
	}

	if ( utility::file::file_extension(resolved_path) == "h5" ) {
#ifdef USEHDF5
		H5FragmentStoreBackend backend(store_path);
		target_store = backend.get_fragment_store(lookup_name);
#else
		utility_exit_with_message("StructureStoreManager::load_fragment_lookup without HDF5 support, unable to load lookup: " + lookup_name + " resolved from: " + store_path);
#endif
	} else {
		BinaryFragmentStoreBackend backend(resolved_path);
		target_store = backend.get_fragment_store(lookup_name);
	}

	FragmentLookupOP lookup( new FragmentLookup(target_store) );
	return lookup;
}

FragmentLookupOP StructureStoreManager::load_fragment_lookup(std::string lookup_name, std::string store_path, std::string group_field,std::string group_type)
{
	FragmentStoreOP target_store(NULL);
	TR.Debug << "loading fragment lookup "+lookup_name+" store path "+store_path+" group_field "+ group_field+ " group type " + group_type;
	std::string resolved_path = resolve_store_path(store_path);

	if ( resolved_path.empty() ) {
		utility_exit_with_message("Unable to resolve specified store: " + store_path);
	}

	if ( utility::file::file_extension(resolved_path) == "h5" ) {
#ifdef USEHDF5
		H5FragmentStoreBackend backend(store_path);
		target_store = backend.get_fragment_store(lookup_name,group_field,group_type);
#else
		utility_exit_with_message("StructureStoreManager::load_fragment_lookup without HDF5 support, unable to load lookup: " + lookup_name + " resolved from: " + store_path);
#endif
	} else {
		BinaryFragmentStoreBackend backend(resolved_path);
		target_store = backend.get_fragment_store(lookup_name);
	}

	FragmentLookupOP lookup = FragmentLookupOP(new FragmentLookup(target_store));
	return lookup;
}

std::map<Size, FragmentStoreOP> StructureStoreManager::group_fragment_store_int(std::string group_field,FragmentStoreOP fullStore){
	//get counts
	std::map <Size,Size> type_ct;
	std::map <Size,Size>::iterator type_ct_iter;
	std::vector<Size> fragsGroupIds = fullStore->int64_groups[group_field];
	Size fragment_length = fullStore->fragment_specification.fragment_length;
	for ( Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
		//second time id seen increment counter 1st time increment to 1
		type_ct_iter = type_ct.find(fragsGroupIds[ii]);
		if ( type_ct_iter != type_ct.end() ) {
			type_ct_iter->second++;
		} else {
			type_ct.insert(std::pair<Size,Size>(fragsGroupIds[ii],1));
		}
	}
	//create correctly sized fragmentStores
	std::map<Size,FragmentStoreOP> typed_frags;
	FragmentSpecification fragment_spec = fullStore->fragment_specification;
	for ( type_ct_iter=type_ct.begin(); type_ct_iter != type_ct.end(); ++type_ct_iter ) {
		FragmentStoreOP fragment_store = FragmentStoreOP(new FragmentStore(fragment_spec, type_ct_iter->second ));
		typed_frags.insert(std::pair<Size,FragmentStoreOP> (type_ct_iter->first, fragment_store));
	}
	//populate fragment stores
	std::map <Size,Size> type_currentCt;
	for ( Size ii =0; ii<fragsGroupIds.size(); ++ii ) {
		Size newFragStartPosition = 0; //position in new fragment array the fragment should start
		Size originalFragStartPosition = ii*fragment_length;
		type_ct_iter = type_currentCt.find(fragsGroupIds[ii]);
		if ( type_ct_iter != type_currentCt.end() ) {
			newFragStartPosition = type_ct_iter->second*fragment_length;
			type_ct_iter->second++;
		} else {
			newFragStartPosition = 0*fragment_length;  //always will be 0 just set so I remember what is actually going on.
			type_currentCt.insert(std::pair<Size,Size>(fragsGroupIds[ii],1));
		}
		for ( Size jj=0; jj<fragment_length; ++jj ) {
			typed_frags[fragsGroupIds[ii]]->fragment_coordinates[newFragStartPosition+jj] =fullStore->fragment_coordinates[originalFragStartPosition+jj];
		}
	}
	return(typed_frags);
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
	if ( utility::file::file_exists(target_path) ) {
		return target_path;
	}

	if ( utility::file::file_exists(basename) ) {
		return basename;
	}

	return "";
}

StructureStoreManager * StructureStoreManager::create_singleton_instance()
{
	return new StructureStoreManager;
}

}
}
