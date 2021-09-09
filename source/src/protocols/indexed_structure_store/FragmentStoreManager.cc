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
#include <utility/file/file_sys_util.hh>

#include <boost/algorithm/string.hpp>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <core/sequence/SSManager.hh>

#include <protocols/indexed_structure_store/FragmentStoreManager.hh>
#include <protocols/indexed_structure_store/FragmentStore.hh>
#include <protocols/indexed_structure_store/FragmentLookup.hh>
#include <protocols/indexed_structure_store/SSHashedFragmentStore.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/StructureStoreManager.hh>

// H5-based backend declarations are guarded by #ifdef USEHDF5
#include <protocols/indexed_structure_store/H5FragmentStoreBackend.hh>
#include <protocols/indexed_structure_store/BinaryFragmentStoreBackend.hh>

#include <utility/vector1.hh>
#include <array>
#include <chrono>
#include <queue>
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
	for ( core::Size ii=1; ii<=fields_to_load.size(); ++ii ) {
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
		std::map <core::Size,Size> type_ct;
		std::map <core::Size,Size>::iterator type_ct_iter;
		std::vector<core::Size> fragsGroupIds = fullStore->int64_groups[group_field];
		core::Size fragment_length = fullStore->fragment_specification.fragment_length;
		for ( core::Size & fragsGroupId : fragsGroupIds ) {
			//second time id seen increment counter 1st time increment to 1
			type_ct_iter = type_ct.find(fragsGroupId);
			if ( type_ct_iter != type_ct.end() ) {
				type_ct_iter->second++;
			} else {
				type_ct.insert(std::pair<core::Size,Size>(fragsGroupId,1));
			}
		}
		//create correctly sized fragmentStores
		std::map<core::Size,FragmentStoreOP> typed_frags;
		FragmentSpecification fragment_spec = fullStore->fragment_specification;
		for ( type_ct_iter=type_ct.begin(); type_ct_iter != type_ct.end(); ++type_ct_iter ) {
			FragmentStoreOP fragment_store = utility::pointer::make_shared< FragmentStore >(fragment_spec, type_ct_iter->second );
			fragment_store->hash_id = type_ct_iter->first;
			typed_frags.insert(std::pair<core::Size,FragmentStoreOP> (type_ct_iter->first, fragment_store));
		}
		//populate fragment stores
		std::map <core::Size,Size> type_currentCt;
		for ( core::Size ii =0; ii<fragsGroupIds.size(); ++ii ) {
			core::Size newFragStartPosition = 0; //position in new fragment array the fragment should start
			core::Size originalFragStartPosition = ii*fragment_length;
			type_ct_iter = type_currentCt.find(fragsGroupIds[ii]);
			if ( type_ct_iter != type_currentCt.end() ) {
				newFragStartPosition = type_ct_iter->second*fragment_length;
				type_ct_iter->second++;
			} else {
				newFragStartPosition = 0*fragment_length;  //always will be 0 just set so I remember what is actually going on.
				type_currentCt.insert(std::pair<core::Size,Size>(fragsGroupIds[ii],1));
			}
			for ( core::Size jj=0; jj<fragment_length; ++jj ) {
				typed_frags[fragsGroupIds[ii]]->fragment_coordinates[newFragStartPosition+jj] =fullStore->fragment_coordinates[originalFragStartPosition+jj];
			}
		}//ends the copying of the fragmentCoordinates
		fullStore->fragment_coordinates.clear();
		vector1<xyzVector<Real> >().swap(fullStore->fragment_coordinates);
		//-------------Copy the real_groups
		for ( auto itr = fullStore->real_groups.begin(); itr !=  fullStore->real_groups.end(); ++itr ) {
			std::string tmp_group_name = itr->first;
			for ( core::Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
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
			for ( core::Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
				if ( typed_frags[fragsGroupIds[ii]]->realVector_groups[tmp_group_name].size() == 0 ) {
					std::vector<std::vector<Real> >real_vector_vector;
					typed_frags[fragsGroupIds[ii]]->realVector_groups.insert(std::pair<std::string,std::vector< std::vector< Real> > > (tmp_group_name,real_vector_vector));
				}
				std::vector<Real> real_vector;
				for ( core::Size jj=0; jj<fullStore->realVector_groups[tmp_group_name][jj].size(); ++jj ) {
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
			for ( core::Size ii=0; ii<fragsGroupIds.size(); ++ii ) {
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
				for ( core::Size ii =0; ii<fragsGroupIds.size(); ++ii ) {
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

//The SSHashedFragmentStore is an instance of the grouped_fragment_store.
SSHashedFragmentStoreOP FragmentStoreManager::SSHashedFragmentStore( std::string store_path,std::string store_format, std::string store_compression)
{
	if ( store_path == "" && SSHashedFragmentStoreMap_.size()>0 ) {
		std::string default_path = SSHashedFragmentStoreMap_.begin()->first;
		if ( store_path != default_path ) { //Do not want to print out if the store_path is "" and the first_path is "" because this person is using flags
			TR << "Auto choosing fragment store :" << store_path<< std::endl;
		}
		return(SSHashedFragmentStoreMap_.begin()->second);
	}
	SSHashedFragmentStoreMap::const_iterator iter(SSHashedFragmentStoreMap_.find(store_path));
	if ( iter != SSHashedFragmentStoreMap_.end() ) {
		return(iter->second);
	} else {
		protocols::indexed_structure_store::SSHashedFragmentStoreOP newHashedFragmentStoreOP( new protocols::indexed_structure_store::SSHashedFragmentStore(store_path,store_format,store_compression) );
		SSHashedFragmentStoreMap_[store_path]=newHashedFragmentStoreOP;
		return(newHashedFragmentStoreOP);
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

std::map<numeric::Size, FragmentStoreOP> FragmentStoreManager::find_previously_loaded_fragmentStore(std::string const & fragment_store_path){
	if ( fragment_store_path.size()<=1 && (grouped_fragment_map_.begin() != grouped_fragment_map_.end()) ) {
		return(grouped_fragment_map_.begin()->second);
	} else {
		std::map<numeric::Size, FragmentStoreOP> empty;
		return empty;
	}

}



std::map<numeric::Size, FragmentStoreOP> FragmentStoreManager::load_and_hash_fragmentStore(std::string const & fragment_store_path, std::string const & fragment_store_compression, numeric::Size const & fragment_size){
	using namespace numeric;
	GroupedFragmentMap::const_iterator iter(grouped_fragment_map_.find(fragment_store_path));
	if ( iter != grouped_fragment_map_.end() ) {
		return(iter->second);
	} else {
		std::map<numeric::Size, FragmentStoreOP> hashed_fragmentStore;
		core::sequence::SSManager SM;
		auto start_time = std::chrono::high_resolution_clock::now();
		StructureStoreOP structure_store = StructureStoreManager::get_instance()->load_structure_store(fragment_store_path);
		//Step 1:  Loop through structure store and figure how many fragments of each dssp type there is
		std::vector<numeric::Size> ss_vector;
		auto post_disk_load_time = std::chrono::high_resolution_clock::now();
		std::map<numeric::Size, std::vector<numeric::Size> > dssp_fragCt = get_dssp_fragCt(structure_store,fragment_store_compression,fragment_size,ss_vector);
		auto post_dssp_time = std::chrono::high_resolution_clock::now();
		//Step 2: For each SS trype
		//create the space
		std::vector<std::string> atoms;
		atoms.push_back("CA");
		FragmentSpecification frag_spec(9,atoms);
		for ( auto& iter: dssp_fragCt ) {
			FragmentStoreOP fragment_store = utility::pointer::make_shared< FragmentStore >(frag_spec, iter.second.size() );
			//for speed predeclare the space for each element in fragment_store
			Size num_fragments = iter.second.size();
			std::vector < std::vector<Real> > phi_v;
			std::vector < std::vector<Real> > psi_v;
			std::vector < std::vector<Real> > omega_v;
			std::vector < std::vector<Real> > cen_v;
			std::vector< std::string> aa_v;
			std::vector< std::string> name_v;
			phi_v.reserve(num_fragments);
			psi_v.reserve(num_fragments);
			omega_v.reserve(num_fragments);
			cen_v.reserve(num_fragments);
			aa_v.reserve(num_fragments);
			name_v.reserve(num_fragments);
			/*phi_v.resize(num_fragments, std::vector<Real>(fragment_size));
			psi_v.resize(num_fragments, std::vector<Real>(fragment_size));
			omega_v.resize(num_fragments, std::vector<Real>(fragment_size));
			cen_v.resize(num_fragments, std::vector<Real>(fragment_size));
			aa_v.resize(num_fragments, std::string(fragment_size,'x'));
			name_v.resize(num_fragments, std::string(5,'x'));
			*/
			/*
			aa_v.resize(num_fragments);
			name_v.resize(num_fragments);
			for(Size ii=0; ii<num_fragments; ++ii){
			aa_v[ii]=std::string(fragment_size,'x');
			name_v[ii]=std::string(5,'x');
			}*/
			fragment_store->realVector_groups.insert(std::pair<std::string,std::vector < std::vector<Real> > > ("phi",phi_v));
			fragment_store->realVector_groups.insert(std::pair<std::string,std::vector < std::vector<Real> > > ("psi",psi_v));
			fragment_store->realVector_groups.insert(std::pair<std::string,std::vector < std::vector<Real> > > ("omega",omega_v));
			fragment_store->realVector_groups.insert(std::pair<std::string,std::vector < std::vector<Real> > > ("cen",cen_v));
			fragment_store->string_groups.insert(std::pair<std::string,std::vector< std::string> > ("aa",aa_v));
			fragment_store->string_groups.insert(std::pair<std::string,std::vector< std::string> > ("name",name_v));
			fragment_store->hash_id = iter.first;
			hashed_fragmentStore.insert(std::pair<core::Size,FragmentStoreOP> (iter.first, fragment_store));

		}
		//load_db----
		//atom coordinates of CA
		//ss_bin
		//phi,psi,omega
		//cen
		//aa
		//fragment threshold distance is initialized & set within the movers
		//name (5 char)
		for ( auto& iter: dssp_fragCt ) {
			numeric::Size frag_index = 0; //The fragment_coordinates are initialized to the proper size. Would this be smart to do for all of my datatypes.
			for ( auto & iter_fragmentStore: iter.second ) {
				uint32_t name_id = structure_store->residue_entries[iter_fragmentStore].structure_id-1;
				std::string name_s= structure_store->structure_entries[name_id].name;
				std::vector<char> aa_v;
				std::vector<Real> phi_v;
				std::vector<Real> psi_v;
				std::vector<Real> omega_v;
				std::vector<Real> cen_v;
				for ( ndarray::Size ii=0; ii < fragment_size; ++ii ) {
					Real phi = structure_store->residue_entries[iter_fragmentStore+ii].bb.phi;
					Real psi = structure_store->residue_entries[iter_fragmentStore+ii].bb.psi;
					Real omega = structure_store->residue_entries[iter_fragmentStore+ii].bb.omega;
					Real cen = structure_store->residue_entries[iter_fragmentStore+ii].cen;
					char aa = structure_store->residue_entries[iter_fragmentStore+ii].sc.aa;
					std::array<float, 3> CA = structure_store->residue_entries[iter_fragmentStore+ii].orient.CA;
					phi_v.push_back(phi);
					psi_v.push_back(psi);
					omega_v.push_back(omega);
					cen_v.push_back(cen);
					aa_v.push_back(aa);
					numeric::xyzVector <numeric::Real> CA_xyz(CA[0],CA[1],CA[2]);
					hashed_fragmentStore[iter.first]->fragment_coordinates[frag_index*fragment_size+ii]=CA_xyz;
				}
				hashed_fragmentStore[iter.first]->realVector_groups["phi"].push_back(phi_v);
				hashed_fragmentStore[iter.first]->realVector_groups["psi"].push_back(psi_v);
				hashed_fragmentStore[iter.first]->realVector_groups["omega"].push_back(omega_v);
				hashed_fragmentStore[iter.first]->realVector_groups["cen"].push_back(cen_v);
				std::string aa_s(aa_v.begin(),aa_v.end());
				hashed_fragmentStore[iter.first]->string_groups["aa"].push_back(aa_s);
				hashed_fragmentStore[iter.first]->string_groups["name"].push_back(name_s);
				frag_index+=1;
			}
		}
		auto post_convert_db_time = std::chrono::high_resolution_clock::now();
		TR << "disk load time: " << std::chrono::duration_cast<std::chrono::seconds>(post_disk_load_time - start_time).count() << " seconds" << std::endl;
		TR << "convert between db time " << std::chrono::duration_cast<std::chrono::seconds>(post_convert_db_time - post_disk_load_time).count() << " seconds" << std::endl;
		TR << "convert dssp time " << std::chrono::duration_cast<std::chrono::seconds>(post_dssp_time - post_disk_load_time).count() << " seconds" << std::endl;
		StructureStoreManager::get_instance()->delete_structure_store(fragment_store_path);
		return(hashed_fragmentStore);
	}
}


std::map<numeric::Size, std::vector<numeric::Size>> FragmentStoreManager::get_dssp_fragCt(StructureStoreOP structure_store,std::string const & fragment_store_compression,numeric::Size const & fragment_size,std::vector<numeric::Size> & ss_vector){
	std::map<numeric::Size, std::vector<numeric::Size>> dssp_fragCt;
	std::queue<char> frag_dssp;
	core::sequence::SSManager SM;
	std::map<numeric::Size, std::vector<numeric::Size>>::iterator it ;
	char ss;
	for ( ndarray::Size i=0; i < structure_store->residue_entries.getSize<0>(); ++i ) {
		bool second_chain_start = false;
		if ( i<structure_store->residue_entries.getSize<0>()-1 && i>1 ) {
			if ( structure_store->residue_entries[i].structure_id != structure_store->residue_entries[i-1].structure_id ) {
				second_chain_start=true;
			}
			if ( structure_store->residue_entries[i-1].chain_ending ) {
				second_chain_start=true;
			}
		}
		if ( !second_chain_start ) {
			ss = structure_store->residue_entries[i].ss;
			//uint32_t name_id = structure_store->residue_entries[i].structure_id-1; //***THE STRUCTURE ID IN THE DATABASE is 1 OFF
			//std::string name= structure_store->structure_entries[name_id].name;
			frag_dssp.push(ss);
			if ( frag_dssp.size()==fragment_size ) {
				std::string frag_string = queue_to_string(frag_dssp);
				bool add_fragment = true;
				if ( fragment_store_compression=="helix_shortLoop" ) {
					core::Size h_ct = std::count(frag_string.begin(), frag_string.end(), 'H'); //helix
					if ( h_ct<4 ) {
						add_fragment = false;
					}
				}
				if ( fragment_store_compression=="sheet_shortLoop" ) {
					core::Size e_ct = std::count(frag_string.begin(), frag_string.end(), 'E'); //sheet
					if ( e_ct<4 ) {
						add_fragment = false;
					}
				}
				if ( add_fragment ) {
					core::Size ss_index = SM.symbolString2index(frag_string);
					ss_vector.push_back(ss_index);
					it = dssp_fragCt.find(ss_index);
					if ( it == dssp_fragCt.end() ) {
						//create vector of positions for new SS type
						std::vector<numeric::Size> tmp_vector;
						tmp_vector.push_back(i-fragment_size+1);
						dssp_fragCt.insert(std::make_pair(ss_index,tmp_vector));
					} else {
						it->second.push_back(i-fragment_size+1);
					}
				}
				//frag_seq.pop();
				frag_dssp.pop();
			}
		} else { //for thse second chain start
			frag_dssp = std::queue<char>(); //empties queue if end of pdb or chain
			ss = structure_store->residue_entries[i].ss;
			frag_dssp.push(ss);
		}
	}
	return(dssp_fragCt);
}

std::string FragmentStoreManager::queue_to_string(std::queue<char> frag_dssp){
	std::string frag_dssp_string ="";
	while ( ! frag_dssp.empty() ) {
		frag_dssp_string +=frag_dssp.front();
		frag_dssp.pop();
	}
	return(frag_dssp_string);
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
