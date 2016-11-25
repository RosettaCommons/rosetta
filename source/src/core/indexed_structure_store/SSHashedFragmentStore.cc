// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @details Parses the VALL by secondary structure after it is read in hdf5 format
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
#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/SSManager.hh>
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


static THREAD_LOCAL basic::Tracer TR( "core.indexed_structure_store.SSHashedFragmentStore" );


namespace core {
namespace indexed_structure_store {
using namespace core;
using utility::vector1;
using namespace std;

SSHashedFragmentStore::SSHashedFragmentStore(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace OptionKeys::indexed_structure_store;
	using namespace core::indexed_structure_store;
	std::string store_path = option[OptionKeys::indexed_structure_store::fragment_store](); //error checking occurs in DB loading
	vector1<string> fields_to_load;
	vector1<string> fields_to_load_types;
	string store_name = option[OptionKeys::indexed_structure_store::store_name]();
	string group_field ="ss_bin";
	fields_to_load.push_back("ss_bin");
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
	if ( option[OptionKeys::indexed_structure_store::exclude_homo].user() ) {
		fields_to_load.push_back("name");
		fields_to_load_types.push_back("five_char_per_residue");
	}
	SSHashedFragmentStore_=StructureStoreManager::get_instance()->load_grouped_fragment_store(group_field,store_name,store_path,fields_to_load,fields_to_load_types);
}

void SSHashedFragmentStore::set_threshold_distance(Real threshold_distance){
	std::map<Size, core::indexed_structure_store::FragmentStoreOP>::iterator fragStoreMap_iter;
	if ( SSHashedFragmentStore_.begin()->second->fragment_threshold_distances[0]>threshold_distance ) { //fragment threshold distance needs to be set to the lowest calling value.
		for ( fragStoreMap_iter = SSHashedFragmentStore_.begin(); fragStoreMap_iter != SSHashedFragmentStore_.end(); fragStoreMap_iter++ ) {
			fragStoreMap_iter->second->add_threshold_distance_allFrag(threshold_distance);
		}
	}
}


void SSHashedFragmentStore::init_SS_stub_HashedFragmentStore(){
	Size fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	vector1<std::string> types;
	types.push_back("HH");
	types.push_back("EE");
	types.push_back("HHH"); //allow flexibility for the end of the helix
	types.push_back("EEE");
	vector1<std::string> loops;
	std::string tmp_loop = "L";
	loops.push_back(tmp_loop);
	vector1<std::string> ss_stubs;
	for ( Size ii=2; ii<=5; ++ii ) {
		tmp_loop+="L";
		loops.push_back(tmp_loop);
	}
	for ( Size ii=1; ii<=types.size(); ++ii ) {
		for ( Size jj=1; jj<=loops.size(); ++jj ) {
			for ( Size kk=1; kk<=types.size(); ++kk ) {
				std::string ss_stub = types[ii]+loops[jj]+types[kk];
				if ( ss_stub.size()<=fragment_length ) {
					ss_stubs.push_back(ss_stub);
				}
			}
		}
	}
	core::sequence::SSManager SM;
	std::map<Size, core::indexed_structure_store::FragmentStoreOP>::iterator fragStoreMap_iter;
	for ( fragStoreMap_iter = SSHashedFragmentStore_.begin(); fragStoreMap_iter != SSHashedFragmentStore_.end(); fragStoreMap_iter++ ) {
		std::string fragStore_ss_string = SM.index2symbolString(fragStoreMap_iter->first,fragment_length);
		for ( Size ii=1; ii<=ss_stubs.size(); ++ii ) {
			if ( fragStore_ss_string.substr(0,ss_stubs[ii].size())==ss_stubs[ii] ) {
				std::vector<Size> residues;
				residues.push_back(0);
				residues.push_back(1);
				residues.push_back(ss_stubs[ii].size()-2);
				residues.push_back(ss_stubs[ii].size()-1);
				fragStoreMap_iter->second->generate_residue_subset_fragment_store(residues);
				if ( SS_stub_HashedFragmentStoreIndex_.find(ss_stubs[ii]) == SS_stub_HashedFragmentStoreIndex_.end() ) {
					vector1<Size> fragStoreIndex_vector;
					SS_stub_HashedFragmentStoreIndex_.insert(std::pair<std::string,vector1<Size> > (ss_stubs[ii],fragStoreIndex_vector));
				}
				SS_stub_HashedFragmentStoreIndex_[ss_stubs[ii]].push_back(fragStoreMap_iter->first);
			}
		}
	}
}


Size SSHashedFragmentStore::get_valid_resid(core::pose::Pose const pose,int resid){
	using namespace core::indexed_structure_store;
	Size fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	if ( resid<1 ) {
		TR << "invalid resid encountered n-term" << resid << std::endl;
		resid = 1;
	}
	if ( resid+fragment_length-1>pose.total_residue() ) {
		TR << "invalid resid encountered c-term" << resid << std::endl;
		resid = (int)pose.total_residue()-fragment_length+1;
	}
	return(Size(resid));
}

set<std::string> SSHashedFragmentStore::potential_valid_ss_strings(std::string frag_ss){
	//dssp sometimes makes errors this corrects some of thems.
	//Note: this was not designed to fix cases where two inaccuracies happened in the same fragment. ie HLLLHLH
	set<std::string> valid_frag_ss;
	set<std::string> tmp_frag_ss;
	set<std::string> union_frag_ss;
	valid_frag_ss.insert(frag_ss);
	set<std::string>::iterator iter;
	//step1 HL -> LL,HH
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		Size found = iter->find("HL");
		while ( found!=std::string::npos ) {
			std::string tmp_string1 = *iter;
			std::string tmp_string2 = *iter;
			tmp_string1[found] ='L';
			tmp_string2[found+1] ='H';
			tmp_frag_ss.insert(tmp_string1);
			tmp_frag_ss.insert(tmp_string2);
			found = iter->find("HL",found+2);
		}
	}
	set_union(valid_frag_ss.begin(),valid_frag_ss.end(),tmp_frag_ss.begin(),tmp_frag_ss.end(),inserter(union_frag_ss,union_frag_ss.begin()));
	valid_frag_ss = union_frag_ss;
	tmp_frag_ss.clear();
	union_frag_ss.clear();
	//step2 EL -> LL,EE
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		set<std::string> union_frag_ss;
		Size found = iter->find("EL");
		while ( found!=std::string::npos ) {
			std::string tmp_string1 = *iter;
			std::string tmp_string2 = *iter;
			tmp_string1[found] ='L';
			tmp_string2[found+1] ='E';
			tmp_frag_ss.insert(tmp_string1);
			tmp_frag_ss.insert(tmp_string2);
			found = iter->find("EL",found+2);
		}
	}
	set_union(valid_frag_ss.begin(),valid_frag_ss.end(),tmp_frag_ss.begin(),tmp_frag_ss.end(),inserter(union_frag_ss,union_frag_ss.begin()));
	valid_frag_ss = union_frag_ss;
	tmp_frag_ss.clear();
	union_frag_ss.clear();
	//step3 LH -> LL,HH
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		set<std::string> union_frag_ss;
		Size found = iter->find("LH");
		while ( found!=std::string::npos ) {
			std::string tmp_string1 = *iter;
			std::string tmp_string2 = *iter;
			tmp_string1[found] ='H';
			tmp_string2[found+1] ='L';
			tmp_frag_ss.insert(tmp_string1);
			tmp_frag_ss.insert(tmp_string2);
			found = iter->find("LH",found+2);
		}
	}
	set_union(valid_frag_ss.begin(),valid_frag_ss.end(),tmp_frag_ss.begin(),tmp_frag_ss.end(),inserter(union_frag_ss,union_frag_ss.begin()));
	valid_frag_ss = union_frag_ss;
	tmp_frag_ss.clear();
	union_frag_ss.clear();
	//step4 LE -> LL,EE
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		set<std::string> union_frag_ss;
		Size found = iter->find("LE");
		while ( found!=std::string::npos ) {
			std::string tmp_string1 = *iter;
			std::string tmp_string2 = *iter;
			tmp_string1[found] ='E';
			tmp_string2[found+1] ='L';
			tmp_frag_ss.insert(tmp_string1);
			tmp_frag_ss.insert(tmp_string2);
			found = iter->find("EL",found+2);
		}
	}
	set_union(valid_frag_ss.begin(),valid_frag_ss.end(),tmp_frag_ss.begin(),tmp_frag_ss.end(),inserter(union_frag_ss,union_frag_ss.begin()));
	valid_frag_ss = union_frag_ss;
	tmp_frag_ss.clear();
	union_frag_ss.clear();
	return(valid_frag_ss);
}


Real SSHashedFragmentStore::max_rmsd_in_region(pose::Pose const pose, vector1<Size> resids){
	Real high_rmsd=0;
	for ( Size ii=1; ii<=resids.size(); ++ii ) {
		Size valid_resid = get_valid_resid(pose,resids[ii]);
		Real tmp_rmsd = lookback_account_for_dssp_inaccuracy(pose,valid_resid,true,0.0);
		if ( tmp_rmsd >high_rmsd ) {
			high_rmsd=tmp_rmsd;
		}
	}
	return(high_rmsd);
}

Real SSHashedFragmentStore::lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid,bool find_closest, Real rms_threshold){
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	Size fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string frag_ss = dssp_string.substr(resid-1,fragment_length);
	return(lookback_account_for_dssp_inaccuracy(pose,resid,frag_ss,find_closest,rms_threshold));
}

Real SSHashedFragmentStore::lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid,std::string frag_ss, bool find_closest, Real rms_threshold){
	set<std::string> valid_frag_ss = potential_valid_ss_strings(frag_ss);
	set<std::string>::iterator iter;
	iter = valid_frag_ss.find("HHHHHHHHH"); //only check all helical if that was requested. This is for time.
	if ( iter!=valid_frag_ss.end() ) {
		valid_frag_ss.erase(iter);
	}
	iter = valid_frag_ss.find(frag_ss);
	if ( iter!=valid_frag_ss.end() ) {
		valid_frag_ss.erase(iter);
	}
	Real low_rmsd = lookback(pose,resid,frag_ss,find_closest); //tries original first.
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		if ( low_rmsd>rms_threshold || find_closest ) {
			Real tmp_rmsd = lookback(pose,resid,*iter,find_closest);
			if ( tmp_rmsd<low_rmsd ) {
				low_rmsd = tmp_rmsd;
			}
		}
	}
	return(low_rmsd);
}

Real SSHashedFragmentStore::lookback_account_for_dssp_inaccuracy(pose::Pose const pose, Size resid, std::string frag_ss, Real & match_rmsd, Size & match_index, Size & match_ss_index){
	using namespace core::indexed_structure_store;
	set<std::string> valid_frag_ss = potential_valid_ss_strings(frag_ss);
	set<std::string>::iterator iter;
	iter = valid_frag_ss.find("HHHHHHHHH"); //only check all helical if that was requested. This is for time.
	if ( iter!=valid_frag_ss.end() ) {
		valid_frag_ss.erase(iter);
	}
	Real low_rmsd = 5;
	core::sequence::SSManager SM;
	Size fragment_length = get_fragment_length();
	//get initial resid
	std::vector< numeric::xyzVector<numeric::Real> > coordinates;
	for ( Size ii = 0;  ii < fragment_length; ++ii ) {
		for ( std::string const & atom_name : get_fragment_store()->fragment_specification.fragment_atoms ) {
			coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
		}
	}
	//find top hit
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		Size ss_index = SM.symbolString2index(*iter);
		if ( SSHashedFragmentStore_.find(ss_index)  != SSHashedFragmentStore_.end() ) {
			FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_.at(ss_index);
			FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
			FragmentLookupResult lookupResult = selected_fragLookupOP->lookup_closest_fragment(&coordinates[0]);
			Real tmp_rmsd = lookupResult.match_rmsd;
			Size tmp_index = lookupResult.match_index;
			if ( tmp_rmsd<low_rmsd ) {
				low_rmsd= tmp_rmsd;
				match_ss_index = ss_index;
				match_index = tmp_index;
				match_rmsd = tmp_rmsd;
			}
		}
	}
	return(low_rmsd);
}


Real SSHashedFragmentStore::lookback(pose::Pose const pose, Size resid){
	using namespace core::indexed_structure_store;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	Size fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string frag_ss = dssp_string.substr(resid-1,fragment_length);
	return(lookback(pose,resid,frag_ss,true));
}


Real SSHashedFragmentStore::lookback(pose::Pose const pose, Size resid,string frag_ss,bool find_closest){
	using namespace core::indexed_structure_store;
	core::sequence::SSManager SM;
	Size fragmentStore_fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	Size ss_index = SM.symbolString2index(frag_ss);
	Real returnRmsd;
	if ( SSHashedFragmentStore_.find(ss_index)  != SSHashedFragmentStore_.end() ) {
		//case where item is found in map;
		FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_.at(ss_index);
		FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
			for ( std::string const & atom_name : selected_fragStoreOP->fragment_specification.fragment_atoms ) {
				coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
			}
		}
		FragmentLookupResult lookupResults;
		if ( find_closest ) {
			lookupResults = selected_fragLookupOP->lookup_closest_fragment(&coordinates[0]);
		} else {
			lookupResults = selected_fragLookupOP->lookup_fragment(&coordinates[0]);
		}
		returnRmsd = lookupResults.match_rmsd;
	} else {
		TR.Debug << "SS not found in map!! given rms 5" << std::endl;
		//case where item is not found in map. ABEGO not found in pdb is bad sign. Give this a 999 rmsd value
		returnRmsd = 5;
	}
	return(returnRmsd);
}

vector <bool> SSHashedFragmentStore::generate_subset_residues_to_compare(Size loop_length,Size fragment_length,bool match_tail){
	vector<bool> tmp;
	tmp.push_back(true);
	tmp.push_back(true);
	for ( Size ii=0; ii<loop_length; ++ii ) {
		tmp.push_back(false);
	}
	tmp.push_back(true);
	tmp.push_back(true);
	for ( Size ii=4+loop_length; ii<fragment_length; ++ii ) {
		tmp.push_back(match_tail);
	}
	return(tmp);
}

void SSHashedFragmentStore::lookback_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, char resTypeBeforeLoop,char resTypeAfterLoop, Size loop_length, Real & match_rmsd, Size & match_index, Size & match_ss_index){
	using namespace core::indexed_structure_store;
	Real low_rmsd = 9999;
	//Size tmp_match_index;
	std::map<std::string, vector1<Size> >::iterator iter;
	Size fragment_length = get_fragment_length();
	for ( iter=SS_stub_HashedFragmentStoreIndex_.begin(); iter != SS_stub_HashedFragmentStoreIndex_.end(); ++iter ) {
		if ( (iter->first.size() == 4+loop_length) && (iter->first[0]==resTypeBeforeLoop) && (iter->first[iter->first.size()-1]==resTypeAfterLoop) ) {
			for ( Size ii=1; ii<=iter->second.size(); ++ii ) {
				FragmentStoreOP stub_fragStoreOP = SSHashedFragmentStore_[iter->second[ii]]->residue_subset_fragment_store;
				FragmentLookupOP stub_fragLookupOP = stub_fragStoreOP->get_fragmentLookup();
				vector<bool> residues_to_compare = generate_subset_residues_to_compare(loop_length,fragment_length,false);
				FragmentLookupResult lookupResults= stub_fragLookupOP->lookup_closest_fragment_subset(&coordinates[0],residues_to_compare);
				Real tmp_rmsd = lookupResults.match_rmsd;
				Size tmp_index = lookupResults.match_index;
				if ( tmp_rmsd<low_rmsd ) {
					low_rmsd= tmp_rmsd;
					match_ss_index = iter->second[ii];
					match_index = tmp_index;
					match_rmsd = tmp_rmsd;
				}
			}
		}
	}
}

void SSHashedFragmentStore::lookback_uncached_stub(std::vector< numeric::xyzVector<numeric::Real> > coordinates, Size stub_match_ss_index, Size loop_length, Real & match_rmsd, Size & match_index){
	//Note: the short_stub_ss is used to minimize the number of abego types to search.
	using namespace core::indexed_structure_store;
	Real low_rmsd = 9999;
	//Size tmp_match_index;
	numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &coordinates.front().x() , coordinates.size());
	numeric::alignment::QCP_Kernel<core::Real> qcp;
	vector1<Real> rot_vector;
	FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_[stub_match_ss_index];
	Size fragment_length = get_fragment_length();
	vector<bool> residues_to_compare = generate_subset_residues_to_compare(loop_length,fragment_length,true);
	for ( Size jj=0; jj<selected_fragStoreOP->num_fragments_; ++jj ) {
		std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates;
		for ( Size kk = 0;  kk < residues_to_compare.size(); ++kk ) {
			if ( residues_to_compare[kk] ) {
				fragCoordinates.push_back(selected_fragStoreOP->fragment_coordinates[jj*9+kk]);
			}
		}
		numeric::alignment::QCP_Kernel<core::Real>::remove_center_of_mass( &fragCoordinates.front().x() , fragCoordinates.size());
		Real tmp_rmsd = qcp.calc_centered_coordinate_rmsd( &coordinates.front().x(), &fragCoordinates.front().x(), fragCoordinates.size(), &rot_vector[1]);
		if ( tmp_rmsd<low_rmsd ) {
			low_rmsd = tmp_rmsd;
			match_index = jj;
			match_rmsd = tmp_rmsd;
		}
	}
}
std::vector< numeric::xyzVector<numeric::Real> > SSHashedFragmentStore::get_fragment_coordinates(Size db_index,Size match_index){
	FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_.at(db_index);
	return(selected_fragStoreOP->get_fragment_coordinates(match_index));
}

// std::vector< numeric::xyzVector<numeric::Real> > SSHashedFragmentStore::lookback_xyz(pose::Pose const pose, Size resid){
//  using namespace core::indexed_structure_store;
//  core::sequence::ABEGOManager AM;
//  typedef numeric::xyzVector<Real> Vec;
//  utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
//  Size fragmentStore_fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
//  std::string fragAbegoStr = "";
//  for ( Size ii=0; ii<fragmentStore_fragment_length; ++ii ) {
//   fragAbegoStr += abegoSeq[resid+ii];
//  }
//  Size base5index = AM.symbolString2base5index(fragAbegoStr);
//  TR.Debug << "fragAbegoStr:" << fragAbegoStr << std::endl;
//  if ( SSHashedFragmentStore_.find(base5index) == SSHashedFragmentStore_.end() ) {
//   utility_exit_with_message("ABEGO not found in map. FAILURE");
//  }
//  //case where item is found in map;
//  FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_.at(base5index);
//  FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
//  std::vector< numeric::xyzVector<numeric::Real> > coordinates;
//  std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates;
//  for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
//   for ( std::string const & atom_name : selected_fragStoreOP->fragment_specification.fragment_atoms ) {
//    coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
//   }
//  }
//  FragmentLookupResult lookupResult = selected_fragLookupOP->lookup_closest_fragment(&coordinates[0]);
//  for ( Size ii = 0;  ii < fragmentStore_fragment_length; ++ii ) {
//   Real xTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].x();
//   Real yTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].y();
//   Real zTmp = selected_fragStoreOP->fragment_coordinates[lookupResult.match_index+ii].z();
//   fragCoordinates.push_back(Vec(xTmp,yTmp,zTmp));
//  }
//  return(fragCoordinates);
// }



// vector<FragmentLookupResult> SSHashedFragmentStore::get_N_fragments(std::string abego_string,Size topNFrags){
//  //get number of fragments
//  FragmentStoreOP selected_fragStoreOP = get_fragment_store(abego_string);
//  vector<Size> chosen_fragments;
//  vector<FragmentLookupResult> lookupResults;
//  //numeric::Size random_index = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
//  for ( Size ii=0; ii<topNFrags && ii<selected_fragStoreOP->num_fragments_; ++ii ) {
//   Size random_frag = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
//   while ( std::find(chosen_fragments.begin(),chosen_fragments.end(),random_frag)!= chosen_fragments.end() )
//     random_frag = numeric::random::rg().random_range(0,selected_fragStoreOP->num_fragments_);
//   chosen_fragments.push_back(random_frag);
//  }
//  for ( Size ii=0; ii<chosen_fragments.size(); ++ii ) {
//   FragmentLookupResult tmpFragInfo;
//   tmpFragInfo.match_index = chosen_fragments[ii];
//   lookupResults.push_back(tmpFragInfo);
//  }
//  return(lookupResults);
// }

// struct less_then_match_rmsd
// {
//  inline bool operator() (const FragmentLookupResult& struct1, const FragmentLookupResult& struct2)
//  {
//   return (struct1.match_rmsd < struct2.match_rmsd);
//  }
// };


// vector<FragmentLookupResult> SSHashedFragmentStore::get_topN_fragments(std::string /*selectionType*/,Size topNFrags, pose::Pose const pose, Size resid,Real rms_threshold,std::string fragAbegoStr){
//  vector<FragmentLookupResult> lookupResults = get_fragments_below_rms(pose,resid,rms_threshold,fragAbegoStr);
//  //sort array based on rms
//  std::sort(lookupResults.begin(), lookupResults.end(),less_then_match_rmsd());
//  vector<FragmentLookupResult> topLookupResults;
//  for ( int ii=0; ii<(int)topNFrags; ++ii ) {
//   topLookupResults.push_back(lookupResults[ii]);
//  }
//  return(topLookupResults);
// }

void  SSHashedFragmentStore::get_hits_below_rms(pose::Pose const pose, Size resid, Real rms_threshold, vector1<vector<Real> > & hits_cen, vector1<Real> & hits_rms, vector1<std::string> & hits_aa){
	//get residue coordinates
	core::sequence::SSManager SM;
	//---------------------------------------------------------------------
	std::vector< numeric::xyzVector<numeric::Real> > coordinates;
	FragmentStoreOP tmp_fragStoreOP = get_fragment_store();
	for ( Size ii = 0;  ii < get_fragment_length(); ++ii ) {
		for ( std::string const & atom_name : tmp_fragStoreOP->fragment_specification.fragment_atoms ) {
			coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
		}
	}
	//get frag dssp
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	std::string dssp_string = dssp.get_dssp_secstruct();
	Size fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string frag_ss = dssp_string.substr(resid-1,fragment_length);
	//generate valid SS types while removing all helical
	set<std::string> valid_frag_ss = potential_valid_ss_strings(frag_ss);
	set<std::string>::iterator iter;
	if ( frag_ss != "HHHHHHHHH" ) {
		iter = valid_frag_ss.find("HHHHHHHHH");
		if ( iter!=valid_frag_ss.end() ) {
			valid_frag_ss.erase(iter);
		}
	}
	//loop through valid SS and collect fragments
	for ( iter=valid_frag_ss.begin(); iter!=valid_frag_ss.end(); ++iter ) {
		Size ss_index = SM.symbolString2index(*iter);
		if ( SSHashedFragmentStore_.find(ss_index)!= SSHashedFragmentStore_.end() ) {
			FragmentStoreOP selected_fragStoreOP = SSHashedFragmentStore_.at(ss_index);
			FragmentLookupOP selected_fragLookupOP = selected_fragStoreOP->get_fragmentLookup();
			vector<FragmentLookupResult> lookupResults = selected_fragLookupOP->lookup_close_fragments(&coordinates[0], rms_threshold);
			for ( auto & lookupResult : lookupResults ) {
				hits_rms.push_back(lookupResult.match_rmsd);
				std::vector<Real> cen_list_frag = selected_fragStoreOP->realVector_groups["cen"][lookupResult.match_index];
				hits_cen.push_back(cen_list_frag);
				std::string tmp_aa = selected_fragStoreOP->string_groups["aa"][lookupResult.match_index];
				hits_aa.push_back(tmp_aa);
			}
		}
	}
}



FragmentStoreOP SSHashedFragmentStore::get_fragment_store(Size db_index){
	if ( SSHashedFragmentStore_.find(db_index) != SSHashedFragmentStore_.end() ) {
		return(SSHashedFragmentStore_.at(db_index));
	}
	return(NULL);
}

FragmentStoreOP SSHashedFragmentStore::get_fragment_store(){
	return(SSHashedFragmentStore_.begin()->second);
}

std::map<Size, core::indexed_structure_store::FragmentStoreOP> SSHashedFragmentStore::get_hashed_fragment_store(){
	return(SSHashedFragmentStore_);
}


Size SSHashedFragmentStore::get_fragment_length(){
	using namespace core::indexed_structure_store;
	Size fragmentStore_fragment_length = SSHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	return(fragmentStore_fragment_length);
}

SSHashedFragmentStore *  SSHashedFragmentStore::create_singleton_instance()
{
	return new SSHashedFragmentStore;
}


}
}

