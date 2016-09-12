// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief A wrapper to allow me to load the hashed fragment store on the datacache.
/// @author TJ Brunette tjbrunette@gmail.com


#include <utility/pointer/ReferenceCount.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/StructureStoreManager.hh>
#include <core/indexed_structure_store/H5FragmentStoreBackend.hh>
#include <core/indexed_structure_store/ABEGOHashedFragmentStore.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <numeric/alignment/QCP_Kernel.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <map>
#include <set>

namespace utility {
#ifdef MULTI_THREADED
template <> std::mutex utility::SingletonBase< core::indexed_structure_store::ABEGOHashedFragmentStore >::singleton_mutex_{};
template <> std::atomic< core::indexed_structure_store::ABEGOHashedFragmentStore * > utility::SingletonBase< core::indexed_structure_store::ABEGOHashedFragmentStore >::instance_( 0 );
#else
template <> core::indexed_structure_store::ABEGOHashedFragmentStore * utility::SingletonBase< core::indexed_structure_store::ABEGOHashedFragmentStore >::instance_( nullptr );
#endif
}

#define foreach BOOST_FOREACH

static THREAD_LOCAL basic::Tracer TR( "core.indexed_structure_store.ABEGOHashedFragmentStore" );


namespace core {
namespace indexed_structure_store {
using namespace core;
using utility::vector1;
using namespace std;

ABEGOHashedFragmentStore::ABEGOHashedFragmentStore(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace OptionKeys::indexed_structure_store;
	using namespace core::indexed_structure_store;
	std::string store_path = option[OptionKeys::indexed_structure_store::fragment_store](); //error checking occurs in DB loading
	vector1<string> fields_to_load;
	vector1<string> fields_to_load_types;
	string store_name = option[OptionKeys::indexed_structure_store::store_name]();
	string group_field ="abego_bin";
	fields_to_load.push_back("abego_bin");
	fields_to_load_types.push_back("int64");
	fields_to_load.push_back("phi");
	fields_to_load_types.push_back("real_per_residue");
	fields_to_load.push_back("psi");
	fields_to_load_types.push_back("real_per_residue");
	fields_to_load.push_back("omega");
	fields_to_load_types.push_back("real_per_residue");
	fields_to_load.push_back("cen");
	fields_to_load_types.push_back("real_per_residue");
	fields_to_load.push_back("aa");
	fields_to_load_types.push_back("char_per_residue");
	fields_to_load.push_back("ss");
	fields_to_load_types.push_back("char_per_residue");
	ABEGOHashedFragmentStore_=StructureStoreManager::get_instance()->load_grouped_fragment_store(group_field,store_name,store_path,fields_to_load,fields_to_load_types);
}


void ABEGOHashedFragmentStore::set_threshold_distance(Real threshold_distance){
	std::map<Size, core::indexed_structure_store::FragmentStoreOP>::iterator fragStoreMap_iter;
	if ( ABEGOHashedFragmentStore_.begin()->second->fragment_threshold_distances[0]>threshold_distance ) { //fragment threshold distance needs to be set to the lowest calling value.
		for ( fragStoreMap_iter = ABEGOHashedFragmentStore_.begin(); fragStoreMap_iter != ABEGOHashedFragmentStore_.end(); fragStoreMap_iter++ ) {
			fragStoreMap_iter->second->add_threshold_distance_allFrag(threshold_distance);
		}
	}
}
void ABEGOHashedFragmentStore::init_ss_stub_to_abego(){
	vector1<std::string> types;
	types.push_back("HH");
	types.push_back("EE");
	vector1<std::string> loops;
	std::string tmp_loop = "L";
	loops.push_back(tmp_loop);
	for ( Size ii=2; ii<=5; ++ii ) {
		tmp_loop+="L";
		loops.push_back(tmp_loop);
	}
	for ( Size ii=1; ii<=types.size(); ++ii ) {
		for ( Size jj=1; jj<=loops.size(); ++jj ) {
			for ( Size kk=1; kk<=types.size(); ++kk ) {
				std::string ss_stub = types[ii]+loops[jj]+types[kk];
				vector1<Size> tmp_vector;
				ss_stub_to_abego_[ss_stub]=tmp_vector;
			}
		}
	}
}

set<std::string> ABEGOHashedFragmentStore::get_ss_stubs_per_fragmentStoreOP(Size base5ABEGOindex, core::indexed_structure_store::FragmentStoreOP fragment_storeOP){
	set<std::string> valid_stubs_present;
	//loop must have valid ABEGO and SS to be considered.
	core::sequence::ABEGOManager AM;
	std::string abego_string = AM.base5index2symbolString(base5ABEGOindex,9);
	if ( abego_string=="AAAAAAAAA" ||abego_string=="BBBBBBBBB" ) {  //Don't consider helix or sheet fragments.
		return valid_stubs_present; //Don't consider AAAAAAAAA fragments with
	}
	for ( Size ii=0; ii<fragment_storeOP->num_fragments_; ++ii ) {//loop through the fragments
		std::string full_string = fragment_storeOP->string_groups["ss"][ii];
		for ( Size jj=5; jj<=9; ++jj ) {//1mer - 5mer plus 2 stubs
			std::string short_string = full_string.substr(0,jj);//1mer loop
			if ( ss_stub_to_abego_.count(short_string) && valid_ss_stub_abego_match(short_string,abego_string) ) {
				valid_stubs_present.insert(short_string);
			}
		}
	}
	return(valid_stubs_present);
}

bool ABEGOHashedFragmentStore::valid_ss_stub_abego_match(std::string ss_stub,std::string abego_string){
	//this checks thath the ss_stub helix residues are actually abego A's and ss_stub sheets are actually abego B's. Else this is an inproperly sized loop.
	bool valid = true;
	for ( Size ii=0; ii<ss_stub.size(); ++ii ) {
		std::string ss = ss_stub.substr(ii,1);
		std::string abego = abego_string.substr(ii,1);
		if ( (ss=="H" && abego!="A") || (ss=="E"&&abego!="B") ) {
			valid=false;
		}
	}
	return valid;
}

void ABEGOHashedFragmentStore::generate_ss_stub_to_abego(){
	init_ss_stub_to_abego();
	std::map<Size, core::indexed_structure_store::FragmentStoreOP>::iterator fragStoreMap_iter;
	for ( fragStoreMap_iter = ABEGOHashedFragmentStore_.begin(); fragStoreMap_iter != ABEGOHashedFragmentStore_.end(); fragStoreMap_iter++ ) {
		set<std::string> valid_ss_stubs = get_ss_stubs_per_fragmentStoreOP(fragStoreMap_iter->first,fragStoreMap_iter->second);
		for ( auto const & valid_ss_stub : valid_ss_stubs ) {
			ss_stub_to_abego_[valid_ss_stub].push_back(fragStoreMap_iter->first);
			if ( fragStoreMap_iter->second->fragmentStore_groups.count(valid_ss_stub)==0 ) { //checking for map existance
				std::vector<Size> residues;
				residues.push_back(0);
				residues.push_back(1);
				residues.push_back(valid_ss_stub.size()-2);
				residues.push_back(valid_ss_stub.size()-1);
				fragStoreMap_iter->second->generate_subset_fragment_store(residues,valid_ss_stub);
			}
		}
	}
	TR << "generated stub database" << std::endl;
}

Real ABEGOHashedFragmentStore::lookback(pose::Pose const pose, Size resid){
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Error not designed to work with symmetric pose");
	}
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string fragAbegoStr = "";
	for ( Size ii=0; ii<fragmentStore_fragment_length; ++ii ) {
		fragAbegoStr += abegoSeq[resid+ii];
	}
	return(lookback(pose,resid,fragAbegoStr));
}

std::vector<bool>  ABEGOHashedFragmentStore::generate_ss_subset_match(FragmentStoreOP fragStoreOP,string stub_ss){
	std::vector<bool> ss_matches;
	for ( Size ii=0; ii<fragStoreOP->num_fragments_; ++ii ) {
		string full_ss = fragStoreOP->string_groups["ss"][ii];
		string short_ss = full_ss.substr(0,stub_ss.size());
		if ( short_ss == stub_ss ) {
			ss_matches.push_back(true);
		} else {
			ss_matches.push_back(false);
		}
	}
	return(ss_matches);
}

void ABEGOHashedFragmentStore::lookback_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string stub_ss, Real & match_rmsd, Size & match_index, Size & match_abego){
	//return rmsd match, abego_type(size) and index
	using namespace core::indexed_structure_store;
	vector1<Size> abego_types = ss_stub_to_abego_[stub_ss];
	Real low_rmsd = 9999;
	//Size tmp_match_index;
	for ( Size ii=1; ii<=abego_types.size(); ii++ ) {
		FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(abego_types[ii]);
		FragmentStoreOP stub_fragStoreOP = selected_fragStoreOP->fragmentStore_groups[stub_ss];
		FragmentLookupOP stub_fragLookupOP = stub_fragStoreOP->get_fragmentLookup();
		vector<bool> ss_subset_match = generate_ss_subset_match(selected_fragStoreOP,stub_ss);
		FragmentLookupResult lookupResults= stub_fragLookupOP->lookup_closest_fragment_subset(&coordinates[0],ss_subset_match);
		Real tmp_rmsd = lookupResults.match_rmsd;
		Size tmp_index = lookupResults.match_index;
		if ( tmp_rmsd<low_rmsd ) {
			low_rmsd= tmp_rmsd;
			match_abego = abego_types[ii];
			match_index = tmp_index;
			match_rmsd = tmp_rmsd;
		}
	}
}

void ABEGOHashedFragmentStore::lookback_uncached_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, std::string short_stub_ss,std::string full_stub_ss, Real & match_rmsd, Size & match_index, Size & match_abego){
	//Note: the short_stub_ss is used to minimize the number of abego types to search.
	using namespace core::indexed_structure_store;
	vector1<Size> abego_types = ss_stub_to_abego_[short_stub_ss];
	Real low_rmsd = 9999;
	//Size tmp_match_index;
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &coordinates.front().x() , coordinates.size());
	numeric::alignment::QCP_Kernel<core::Real> qcp;
	vector1<Real> rot_vector;
	for ( Size ii=1; ii<=abego_types.size(); ii++ ) {
		FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(abego_types[ii]);
		std::vector<bool> frags_to_check = generate_ss_subset_match(selected_fragStoreOP,short_stub_ss);
		for ( Size jj=0; jj<selected_fragStoreOP->num_fragments_; ++jj ) {
			if ( frags_to_check[jj] ) { //correct secondary structure of stub. Quick filter to limit computation
				std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates;
				for ( Size kk = 0;  kk < full_stub_ss.size(); ++kk ) {
					if ( full_stub_ss.substr(kk,1)!="L" ) {
						fragCoordinates.push_back(selected_fragStoreOP->fragment_coordinates[jj*9+kk]);
					}
				}
				numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &fragCoordinates.front().x() , fragCoordinates.size());
				Real tmp_rmsd = qcp.calc_centered_coordinate_rmsd( &coordinates.front().x(), &fragCoordinates.front().x(), fragCoordinates.size(), &rot_vector[1]);
				if ( tmp_rmsd<low_rmsd ) {
					low_rmsd = tmp_rmsd;
					match_abego = abego_types[ii];
					match_index = jj;
					match_rmsd = tmp_rmsd;
				}
			}
		}
	}
}


Real ABEGOHashedFragmentStore::lookback(pose::Pose const pose, Size resid,string fragAbegoStr){
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	Size base5index = AM.symbolString2base5index(fragAbegoStr);
	Real returnRmsd;
	if ( ABEGOHashedFragmentStore_.find(base5index)  != ABEGOHashedFragmentStore_.end() ) {
		//case where item is found in map;
		FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(base5index);
		FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
			BOOST_FOREACH ( std::string atom_name, selected_fragStoreOP->fragment_specification.fragment_atoms ) {
				coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
			}
		}
		FragmentLookupResult lookupResults= selected_fragLookupOP->lookup_fragment(&coordinates[0]);
		returnRmsd = lookupResults.match_rmsd;
	} else {
		TR.Debug << "ABEGO not found in map!! given rms 9999" << std::endl;
		//case where item is not found in map. ABEGO not found in pdb is bad sign. Give this a 999 rmsd value
		returnRmsd = 999;
	}
	return(returnRmsd);
}

std::vector< numeric::xyzVector<numeric::Real> > ABEGOHashedFragmentStore::lookback_xyz(pose::Pose const pose, Size resid){
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Error not designed to work with symmetric pose");
	}
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string fragAbegoStr = "";
	for ( Size ii=0; ii<fragmentStore_fragment_length; ++ii ) {
		fragAbegoStr += abegoSeq[resid+ii];
	}
	Size base5index = AM.symbolString2base5index(fragAbegoStr);
	TR.Debug << "fragAbegoStr:" << fragAbegoStr << std::endl;
	if ( ABEGOHashedFragmentStore_.find(base5index) == ABEGOHashedFragmentStore_.end() ) {
		utility_exit_with_message("ABEGO not found in map. FAILURE");
	}
	//case where item is found in map;
	FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(base5index);
	FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
	std::vector< numeric::xyzVector<numeric::Real> > coordinates;
	std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates;
	for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
		BOOST_FOREACH ( std::string atom_name, selected_fragStoreOP->fragment_specification.fragment_atoms ) {
			coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
		}
	}
	FragmentLookupResult lookupResult = selected_fragLookupOP->lookup_closest_fragment(&coordinates[0]);
	for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
		Real xTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].x();
		Real yTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].y();
		Real zTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].z();
		fragCoordinates.emplace_back(xTmp,yTmp,zTmp);
	}
	return(fragCoordinates);
}

std::vector< numeric::xyzVector<numeric::Real> > ABEGOHashedFragmentStore::get_fragment_coordinates(Size match_abego,Size match_index){
	FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(match_abego);
	return(selected_fragStoreOP->get_fragment_coordinates(match_index));
}


vector<FragmentLookupResult> ABEGOHashedFragmentStore::get_N_fragments(std::string abego_string,Size topNFrags){
	//get number of fragments
	FragmentStoreOP selected_fragStoreOP = get_fragment_store(abego_string);
	vector<Size> chosen_fragments;
	vector<FragmentLookupResult> lookupResults;
	//numeric::Size random_index = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
	for ( Size ii=0; ii<topNFrags && ii<selected_fragStoreOP->num_fragments_; ++ii ) {
		Size random_frag = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
		while ( std::find(chosen_fragments.begin(),chosen_fragments.end(),random_frag)!= chosen_fragments.end() )
				random_frag = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
		chosen_fragments.push_back(random_frag);
	}
	for ( unsigned long chosen_fragment : chosen_fragments ) {
		FragmentLookupResult tmpFragInfo;
		tmpFragInfo.match_index = chosen_fragment;
		lookupResults.push_back(tmpFragInfo);
	}
	return(lookupResults);
}

struct less_then_match_rmsd
{
	inline bool operator() (const FragmentLookupResult& struct1, const FragmentLookupResult& struct2)
	{
		return (struct1.match_rmsd < struct2.match_rmsd);
	}
};


vector<FragmentLookupResult> ABEGOHashedFragmentStore::get_topN_fragments(std::string /*selectionType*/,Size topNFrags, pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr){
	vector<FragmentLookupResult> lookupResults = get_fragments_below_rms(pose,resid,rms_threshold,fragAbegoStr);
	//sort array based on rms
	std::sort(lookupResults.begin(), lookupResults.end(),less_then_match_rmsd());
	vector<FragmentLookupResult> topLookupResults;
	for ( int ii=0; ii<(int)topNFrags; ++ii ) {
		topLookupResults.push_back(lookupResults[ii]);
	}
	return(topLookupResults);
}



vector<FragmentLookupResult> ABEGOHashedFragmentStore::get_fragments_below_rms(pose::Pose const pose, Size resid, Real rms_threshold){
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Error not designed to work with symmetric pose");
	}
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string fragAbegoStr = "";
	for ( Size ii=0; ii<fragmentStore_fragment_length; ++ii ) {
		fragAbegoStr += abegoSeq[resid+ii];
	}
	return(get_fragments_below_rms(pose, resid,rms_threshold,fragAbegoStr));
}

vector<FragmentLookupResult> ABEGOHashedFragmentStore::get_fragments_below_rms(pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr){
	using namespace core::indexed_structure_store;
	vector<FragmentLookupResult> lookupResults;
	vector1<string> top_hits_aa;
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	core::sequence::ABEGOManager AM;
	Size base5index = AM.symbolString2base5index(fragAbegoStr);
	if ( ABEGOHashedFragmentStore_.find(base5index)  != ABEGOHashedFragmentStore_.end() ) {
		//case where item is found in map;
		FragmentStoreOP selected_fragStoreOP = ABEGOHashedFragmentStore_.at(base5index);
		FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
			BOOST_FOREACH ( std::string atom_name, selected_fragStoreOP->fragment_specification.fragment_atoms ) {
				coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
			}
		}
		lookupResults= selected_fragLookupOP->lookup_close_fragments(&coordinates[0], rms_threshold);
	} else {
		TR.Debug << "ABEGO not found in map!! No sequences used" << std::endl;
		//case where item is not found in map. ABEGO not found in pdb is bad sign. Give this a 999 rmsd value
	}
	return(lookupResults);
}

std::string ABEGOHashedFragmentStore::get_abego_string(pose::Pose const pose, Size resid){
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Error not designed to work with symmetric pose");
	}
	std::string fragAbegoStr = "";
	for ( Size ii=0; ii<fragmentStore_fragment_length; ++ii ) {
		fragAbegoStr += abegoSeq[resid+ii];
	}
	return(fragAbegoStr);
}


FragmentStoreOP ABEGOHashedFragmentStore::get_fragment_store(std::string fragment_abego){
	core::sequence::ABEGOManager AM;
	Size base5index = AM.symbolString2base5index(fragment_abego);
	if ( ABEGOHashedFragmentStore_.find(base5index) != ABEGOHashedFragmentStore_.end() ) {
		return(ABEGOHashedFragmentStore_.at(base5index));
	}
	return(nullptr);
}

FragmentStoreOP ABEGOHashedFragmentStore::get_fragment_store(Size base5index){
	if ( ABEGOHashedFragmentStore_.find(base5index) != ABEGOHashedFragmentStore_.end() ) {
		return(ABEGOHashedFragmentStore_.at(base5index));
	}
	return(nullptr);
}

FragmentStoreOP ABEGOHashedFragmentStore::get_fragment_store(){
	core::sequence::ABEGOManager AM;
	return(ABEGOHashedFragmentStore_.at(0));
}

Size ABEGOHashedFragmentStore::get_fragment_length(){
	using namespace core::indexed_structure_store;
	Size fragmentStore_fragment_length = ABEGOHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	return(fragmentStore_fragment_length);
}

ABEGOHashedFragmentStore *  ABEGOHashedFragmentStore::create_singleton_instance()
{
	return new ABEGOHashedFragmentStore;
}


}
}

