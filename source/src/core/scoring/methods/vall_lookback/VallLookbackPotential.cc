// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VallLookbackPotenial .cc
/// @brief  threshold based distance metric that measures from any fragment to the structure in the vall
/// @author TJ Brunette (tjbrunette@gmail.com)
///
///
#include <core/scoring/methods/vall_lookback/VallLookbackPotential.hh>
#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/StructureStoreManager.hh>
#include <core/indexed_structure_store/H5FragmentStoreBackend.hh>
#include <core/indexed_structure_store/FragmentLookup.hh>

#include <core/scoring/methods/vall_lookback/VallLookbackData.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


#include <core/pose/Pose.hh>

#include <core/types.hh>
#include <core/sequence/ABEGOManager.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

#define foreach BOOST_FOREACH

#include <map>

static basic::Tracer TR( "core.scoring.vall_lookback.VallLookbackPotential" );

namespace core {
namespace scoring {
namespace methods {

VallLookbackPotential::VallLookbackPotential(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace OptionKeys::indexed_structure_store;
	using namespace core::indexed_structure_store;
	std::string db = option[fragment_store]();
	if (!option[fragment_store].user())
		utility_exit_with_message("Must have option indexed_structure_store::fragment_store defined to use the vall lookback score function");
	TR <<  "loading db:" << db << std::endl;
	FragmentStoreOP target_store = NULL;
	#ifdef USEHDF5
		H5FragmentStoreBackend backend(db);
		target_store = backend.get_fragment_store("9_mer","abego_bin","int64");
	#else
		utility_exit_with_message("StructureStoreManager::load_fragment_store without HDF5 support, unable to load store");
	#endif
	std::map<Size, core::indexed_structure_store::FragmentStoreOP>::iterator fragStoreMap_iter;
	//add distance to all fragments
	Real thresholdDistance = option[fragment_threshold_distance]();
	TR << "setting fragment threshold distance to" << thresholdDistance;
	abegoHashedFragmentStore_=StructureStoreManager::get_instance()->group_fragment_store_int("abego_bin",target_store);
	for(fragStoreMap_iter = abegoHashedFragmentStore_.begin(); fragStoreMap_iter != abegoHashedFragmentStore_.end(); fragStoreMap_iter++ ){
		fragStoreMap_iter->second->add_threshold_distance_allFrag(thresholdDistance);
	}
}

Real VallLookbackPotential::lookback(pose::Pose & pose, Size resid) const{
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;
	runtime_assert( resid<=pose.total_residue());
	if(!pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA)){
		VallLookbackDataOP history = VallLookbackDataOP(new VallLookbackData(pose));
		pose.data().set( CacheableDataType::VALL_LOOKBACK_DATA,history);
	}
	VallLookbackDataOP history_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >( pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));

	if(!history_->get_res_changed(resid)){
		return(history_->get_rmsd(resid));
	}
	else{
		return(lookback_db(pose,resid));
	}
}

Real VallLookbackPotential::lookback(const pose::Pose & pose, Size resid) const{
	//const by removing the creation of the datacache. It'll save resources if it exists. But don't buffer the data if not used.
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;
	runtime_assert( resid<=pose.total_residue());
	return(lookback_db(pose,resid));
}


Real VallLookbackPotential::lookback_db(pose::Pose & pose, Size resid) const{
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;
	using namespace core::indexed_structure_store;
	VallLookbackDataOP history_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >( pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	Size fragmentStore_fragment_length = abegoHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string fragAbegoStr = "";
	for(Size ii=0; ii<fragmentStore_fragment_length; ++ii)
		fragAbegoStr += abegoSeq[resid+ii];
	if(fragAbegoStr == "AAAAAAAAA") {
		history_->set_rmsd(resid,0);
		history_->set_res_changed(resid,false);
		return 0;
	}
	Size base5index = AM.symbolString2base5index(fragAbegoStr);
	Real returnRmsd;
	if(abegoHashedFragmentStore_.find(base5index)  != abegoHashedFragmentStore_.end()){
		//case where item is found in map;
		FragmentStoreOP selected_fragStoreOP = abegoHashedFragmentStore_.at(base5index);
		FragmentLookupOP selected_fragLookupOP = FragmentLookupOP(new FragmentLookup(selected_fragStoreOP));
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for (Size ii = 0;  ii < fragmentStore_fragment_length; ++ii){
			BOOST_FOREACH(std::string atom_name, selected_fragStoreOP->fragment_specification.fragment_atoms){
				coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
			}
		}
		FragmentLookupResult lookupResults= selected_fragLookupOP->lookup_fragment(&coordinates[0]);
		history_->set_rmsd(resid, lookupResults.match_rmsd);
		history_->set_res_changed(resid, false);
		returnRmsd = lookupResults.match_rmsd;
	}
	else{
		TR.Debug << "ABEGO not found in map!! given rms 9999" << std::endl;
		//case where item is not found in map. ABEGO not found in pdb is bad sign. Give this a 999 rmsd value
		returnRmsd = 999;
		history_->set_rmsd(resid, returnRmsd);
		history_->set_res_changed(resid, false);
	}
	return(returnRmsd);
}


Real VallLookbackPotential::lookback_db(const pose::Pose & pose, Size resid) const{
	//maintains const by not adding any history to the the datacache. This should not be the version most of the code uses.
	using namespace core::indexed_structure_store;
	core::sequence::ABEGOManager AM;
	utility::vector1< std::string > abegoSeq = AM.get_symbols( pose,1 );//1 stands for class of ABEGO strings
	Size fragmentStore_fragment_length = abegoHashedFragmentStore_.begin()->second->fragment_specification.fragment_length;
	std::string fragAbegoStr = "";
	for(Size ii=0; ii<fragmentStore_fragment_length; ++ii)
		fragAbegoStr += abegoSeq[resid+ii];
	if(fragAbegoStr == "AAAAAAAAA") { //don't bother looking for helical rmsd. All are very close by rms
		return 0;
	}
	Size base5index = AM.symbolString2base5index(fragAbegoStr);
	if(abegoHashedFragmentStore_.find(base5index)  != abegoHashedFragmentStore_.end()){
		FragmentStoreOP selected_fragStoreOP = abegoHashedFragmentStore_.at(base5index);
		FragmentLookupOP selected_fragLookupOP = FragmentLookupOP(new FragmentLookup(selected_fragStoreOP));
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for (Size ii = 0;  ii < fragmentStore_fragment_length; ++ii){
			BOOST_FOREACH(std::string atom_name, selected_fragStoreOP->fragment_specification.fragment_atoms){
				coordinates.push_back(pose.residue(resid+ii).xyz(atom_name));
			}
		}
		FragmentLookupResult lookupResults= selected_fragLookupOP->lookup_fragment(&coordinates[0]);
		return(lookupResults.match_rmsd);
	}
	else{//case where item is not found in map. ABEGO not found in pdb is bad sign. Give this a 999 rmsd value
		return(999);
	}
}


Real VallLookbackPotential::lookbackRegion(pose::Pose & pose, Size startRes, Size endRes) const{
	using namespace core::pose::datacache;
	Real maxRmsd = 0;
	for (Size ii = startRes;  ii <= endRes; ++ii){
		Real tmpRmsd = lookback(pose,ii);
		if(tmpRmsd > maxRmsd)
			maxRmsd = tmpRmsd;
	}
	return(maxRmsd);
}

Real VallLookbackPotential::lookbackRegion(const pose::Pose & pose, Size startRes, Size endRes) const{
	using namespace core::pose::datacache;
	Real maxRmsd = 0;
	for (Size ii = startRes;  ii <= endRes; ++ii){
		Real tmpRmsd = lookback(pose,ii);
		if(tmpRmsd > maxRmsd)
			maxRmsd = tmpRmsd;
	}
	return(maxRmsd);
}


Real VallLookbackPotential::lookback(pose::Pose & pose) const{
	return(lookbackRegion(pose,1,pose.total_residue()-fragLengthInDB()+1));
}

Real VallLookbackPotential::lookback(const pose::Pose & pose) const{
	return(lookbackRegion(pose,1,pose.total_residue()-fragLengthInDB()+1));
}

Size VallLookbackPotential::fragLengthInDB() const {
	return(abegoHashedFragmentStore_.begin()->second->fragment_specification.fragment_length);
}

}//end methods
}//end scoring
}//end core
