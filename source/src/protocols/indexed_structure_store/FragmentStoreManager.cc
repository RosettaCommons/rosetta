// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @file
/// @author Alex Ford <fordas@uw.edu>

#include <utility/exit.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <protocols/indexed_structure_store/FragmentStoreManager.hh>
#include <protocols/indexed_structure_store/FragmentStore.hh>
#include <protocols/indexed_structure_store/FragmentLookup.hh>

// H5-based backend declarations are guarded by #ifdef USEHDF5
#include <protocols/indexed_structure_store/H5FragmentStoreBackend.hh>
#include <protocols/indexed_structure_store/BinaryFragmentStoreBackend.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace indexed_structure_store {

static basic::Tracer TR( "core.indexed_structure_store.FragmentStoreManager" );

using namespace basic::options;
using utility::vector1;

FragmentStoreManager::FragmentStoreManager()  : store_providers()
{
	typedef std::tuple<core::SSize, std::string> Key;
#ifdef USEHDF5
		store_providers[Key(0, "hdf5")] = FragmentStoreProviderOP(new H5FragmentStoreBackend());
#endif
	store_providers[Key(10, "binary")] = utility::pointer::make_shared< BinaryFragmentStoreBackend >();
}

FragmentLookupOP FragmentStoreManager::load_fragment_lookup(std::string lookup_name)
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

//considered caching at this level but I would get double storage
FragmentLookupOP FragmentStoreManager::load_fragment_lookup(std::string lookup_name, std::string store_path)
{
	FragmentStoreOP target_store(nullptr);

	std::string resolved_path = resolve_store_path(store_path);

	if ( resolved_path.empty() ) {
		utility_exit_with_message("Unable to resolve specified store: " + store_path);
	}

	for ( auto & provider : store_providers ) {
		auto & prio = std::get<0>(provider.first);
		auto & name = std::get<1>(provider.first);
		auto backend = provider.second;
		TR.Debug << "Checking backend: " << name << " prio: " << prio << std::endl;
		backend->set_target_filename(store_path);

		FragmentStoreOP target_store = backend->get_fragment_store(lookup_name);
		FragmentLookupOP lookup( new FragmentLookup(target_store) );
		return lookup;
	}

	if ( utility::file::file_extension(store_path) == "h5" ) {
#ifndef USEHDF5
		utility_exit_with_message("StructureStoreManager::load_structure_store without HDF5 support, unable to load: " + store_path);
#endif
	}
	utility_exit_with_message("Unable to load specified store: " + store_path);
}

//considered caching at this level but I would get double storage
FragmentStoreOP FragmentStoreManager::load_fragment_store(std::string lookup_name, std::string store_path, vector1<std::string> fields_to_load, vector1<std::string> fields_to_load_types){
	std::string resolved_path = resolve_store_path(store_path);
	if ( resolved_path.empty() ) {
		utility_exit_with_message("Unable to resolve specified store: " + store_path);
		TR << fields_to_load.size() << "," << fields_to_load_types.size() << std::endl; //For the testing server. Line should never be accessed.
	}
	for ( auto & provider : store_providers ) {
		auto & prio = std::get<0>(provider.first);
		auto & name = std::get<1>(provider.first);
		auto backend = provider.second;
		TR.Debug << "Checking backend: " << name << " prio: " << prio << std::endl;
		backend->set_target_filename(store_path);

		FragmentStoreOP target_store = backend->get_fragment_store(lookup_name);
		for ( numeric::Size ii=1; ii<=fields_to_load.size(); ++ii ) {
			backend->append_to_fragment_store(target_store,lookup_name,fields_to_load[ii],fields_to_load_types[ii]);
		}
		TR << "database loaded!" << std::endl;
		return target_store;
	}
	utility_exit_with_message("No valid FragmentStoreProvider found for file " + store_path);
}


std::map<numeric::Size, FragmentStoreOP> FragmentStoreManager::load_grouped_fragment_store(std::string group_field, std::string lookup_name, std::string store_path, vector1<std::string> fields_to_load, vector1<std::string> fields_to_load_types){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace OptionKeys::indexed_structure_store;
	using namespace numeric;
	std::string if_cached_name;
	if_cached_name = group_field + lookup_name;
	for ( Size ii=1; ii<=fields_to_load.size(); ++ii ) {
		if_cached_name += fields_to_load[ii];
	}
	GroupedFragmentMap::const_iterator iter(grouped_fragment_map_.find(if_cached_name));
	if ( iter != grouped_fragment_map_.end() ) {
		return(iter->second);
	} else {
		FragmentStoreOP fullStore = load_fragment_store(lookup_name,store_path,fields_to_load,fields_to_load_types);
		if ( option[OptionKeys::indexed_structure_store::exclude_homo].user() ) {
			fullStore->delete_homologs();
		}
		// fullStore->fragment_threshold_distances.clear();
		// fullStore->fragment_coordinates.clear();
		// fullStore->int64_groups.clear();
		// fullStore->real_groups.clear();
		// fullStore->realVector_groups.clear();
		// fullStore->string_groups.clear();
		// std::cout <<"Pause hereC-after delete to watch memory memory" << std::endl;
		// usleep(10000000);
		// std::cout << "end pauseC-after delete" << std::endl;
		std::map <Size,Size> type_ct;
		std::map <Size,Size>::iterator type_ct_iter;
		std::vector<Size> fragsGroupIds = fullStore->int64_groups[group_field];
		Size fragment_length = fullStore->fragment_specification.fragment_length;
		for ( core::Size & fragsGroupId : fragsGroupIds ) {
			//second time id seen increment counter 1st time increment to 1
			type_ct_iter = type_ct.find(fragsGroupId);
			if ( type_ct_iter != type_ct.end() ) {
				type_ct_iter->second++;
			} else {
				type_ct.insert(std::pair<Size,Size>(fragsGroupId,1));
			}
		}
		//create correctly sized fragmentStores
		std::map<Size,FragmentStoreOP> typed_frags;
		FragmentSpecification fragment_spec = fullStore->fragment_specification;
		for ( type_ct_iter=type_ct.begin(); type_ct_iter != type_ct.end(); ++type_ct_iter ) {
			FragmentStoreOP fragment_store = utility::pointer::make_shared< FragmentStore >(fragment_spec, type_ct_iter->second );
			fragment_store->hash_id = type_ct_iter->first;
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
		}//ends the copying of the fragmentCoordinates
		fullStore->fragment_coordinates.clear();
		vector1<xyzVector<Real> >().swap(fullStore->fragment_coordinates);
		//-------------Copy the real_groups
		for ( auto itr = fullStore->real_groups.begin(); itr !=  fullStore->real_groups.end(); ++itr ) {
			std::string tmp_group_name = itr->first;
			for ( Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
				if ( typed_frags[fragsGroupIds[ii]]->real_groups[tmp_group_name].size() == 0 ) {
					std::vector<Real> real_vector;
					typed_frags[fragsGroupIds[ii]]->real_groups.insert(std::pair<std::string,std::vector<Real> > (tmp_group_name,real_vector));
				}
				typed_frags[fragsGroupIds[ii]]->real_groups[tmp_group_name].push_back(fullStore->real_groups[tmp_group_name][ii]);
			}
		}
		fullStore->real_groups.clear();
		std::map<std::string, std::vector<Real> >().swap(fullStore->real_groups);
		//-------------Copy the realVector_groups
		for ( auto itr = fullStore->realVector_groups.begin(); itr !=  fullStore->realVector_groups.end(); ++itr ) {
			std::string tmp_group_name = itr->first;
			for ( Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
				if ( typed_frags[fragsGroupIds[ii]]->realVector_groups[tmp_group_name].size() == 0 ) {
					std::vector<std::vector<Real> >real_vector_vector;
					typed_frags[fragsGroupIds[ii]]->realVector_groups.insert(std::pair<std::string,std::vector< std::vector< Real> > > (tmp_group_name,real_vector_vector));
				}
				std::vector<Real> real_vector;
				for ( Size jj=0; jj<fullStore->realVector_groups[tmp_group_name][jj].size(); ++jj ) {
					real_vector.push_back(fullStore->realVector_groups[tmp_group_name][ii][jj]);
				}
				typed_frags[fragsGroupIds[ii]]->realVector_groups[tmp_group_name].push_back(real_vector);
			}
		}
		fullStore->realVector_groups.clear();
		std::map<std::string, std::vector<std::vector<Real> > >().swap(fullStore->realVector_groups);
		//-------------Copy the string_groups
		for ( auto itr = fullStore->string_groups.begin(); itr !=  fullStore->string_groups.end(); ++itr ) {
			std::string tmp_group_name = itr->first;
			for ( Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
				if ( typed_frags[fragsGroupIds[ii]]->string_groups[tmp_group_name].size() == 0 ) {
					std::vector<std::string> string_vector;
					typed_frags[fragsGroupIds[ii]]->string_groups.insert(std::pair<std::string,std::vector< std::string > > (tmp_group_name,string_vector));
				}
				typed_frags[fragsGroupIds[ii]]->string_groups[tmp_group_name].push_back(fullStore->string_groups[tmp_group_name][ii]);
			}
		}
		fullStore->string_groups.clear();
		std::map<std::string, std::vector<std::string> >().swap(fullStore->string_groups);
		//-------------Copy & Erase the int64_groups while ignoring group_field
		for ( auto itr = fullStore->int64_groups.begin(); itr !=  fullStore->int64_groups.end(); ++itr ) {
			std::string tmp_group_name = itr->first;
			if ( tmp_group_name != group_field ) {
				for ( Size ii =0; ii<fragsGroupIds.size(); ++ii ) {
					if ( typed_frags[fragsGroupIds[ii]]->int64_groups[tmp_group_name].size() == 0 ) {
						std::vector<numeric::Size> int64_vector;
						typed_frags[fragsGroupIds[ii]]->int64_groups.insert(std::pair<std::string,std::vector< numeric::Size > > (tmp_group_name,int64_vector));
					}
					typed_frags[fragsGroupIds[ii]]->int64_groups[tmp_group_name].push_back(fullStore->int64_groups[tmp_group_name][ii]);
				}
			}
		}
		fullStore->int64_groups.clear();
		std::map<std::string, std::vector<numeric::Size> >().swap(fullStore->int64_groups);
		fullStore->fragment_threshold_distances.clear();
		std::vector<numeric::Real>().swap(fullStore->fragment_threshold_distances);
		grouped_fragment_map_[if_cached_name]=typed_frags;
		return(typed_frags);
	}
}


std::string FragmentStoreManager::resolve_store_path(std::string target_path)
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


void FragmentStoreManager::register_store_provider(core::SSize priority, std::string name, FragmentStoreProviderOP backend) {

	std::lock_guard<std::mutex> cache_lock(cache_mutex);

	for ( auto & provider : store_providers ) {
		auto & existing_name = std::get<1>(provider.first);
		if ( name == existing_name ) {
			utility_exit_with_message("Attempting to register duplicate provider name in FragmentStoreManager.");
		}

		store_providers[std::make_tuple(priority, name)] = backend;
	}
}


}

}
