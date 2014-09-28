// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterSet.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/CDRClusterSet.hh>

#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>


namespace protocols {
namespace antibody {
namespace clusters {
	using namespace protocols::antibody;
	using namespace basic::datacache;
	
CDRClusterSet::CDRClusterSet(AntibodyInfo * ab_info){
	ab_info_ = ab_info;
	cluster_matcher_ = CDRClusterMatcherOP( new CDRClusterMatcher() );
	clear();
}

CDRClusterSet::~CDRClusterSet(){}

void
CDRClusterSet::identify_and_set_cdr_cluster(core::pose::Pose const & pose, CDRNameEnum cdr){
	clear(cdr);
	
	core::Size start = ab_info_->get_CDR_start(cdr, pose, North);
	core::Size end = ab_info_->get_CDR_end(cdr, pose, North);
	
	CDRClusterOP cluster = cluster_matcher_->get_cdr_cluster(pose, cdr, start, end);
	clusters_[ cdr ] = cluster;
}

void
CDRClusterSet::clear(){
	clusters_.resize(6, NULL);
}

void
CDRClusterSet::clear(CDRNameEnum cdr){
	clusters_[cdr] = NULL;
}

CDRClusterCOP
CDRClusterSet::get_cluster_data(CDRNameEnum cdr) const{
	return clusters_[cdr];
}

void
CDRClusterSet::set_cluster_data(CDRNameEnum cdr, CDRClusterCOP cluster) {
	clusters_[cdr] = cluster->clone();
}

CDRClusterEnum
CDRClusterSet::get_cluster(CDRNameEnum cdr) const {
	if (clusters_[cdr]){
		return clusters_[cdr]->cluster();
	}
	else {
		return NA;
	}
}

bool
CDRClusterSet::empty(CDRNameEnum cdr) const {
	if (clusters_[cdr]){
		return false;
	}
	else {
		return true;
	}
}

bool
CDRClusterSet::empty() const {
	
	for (core::SSize i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (! empty(cdr)){
			return false;
		}
	}
	return true;
}

BasicCDRClusterSetOP
CDRClusterSet::get_cacheable_cluster_data() const {
	//Make sure we are cop
	return BasicCDRClusterSetOP( new BasicCDRClusterSet(clusters_) );
}

void
CDRClusterSet::set_cacheable_cluster_data_to_pose(core::pose::Pose& pose) const {
	using basic::datacache::DataCache_CacheableData;
	pose.data().set(core::pose::datacache::CacheableDataType::CDR_CLUSTER_INFO, DataCache_CacheableData::DataOP( new BasicCDRClusterSet(clusters_) ));
}

////////////////////////////////////////////////////////////////////////////////
//Note - this is not a base class of CDRClusterSet as we do not want CDRClusterSet to be casheable data due to AntibodyInfoAP.
BasicCDRClusterSet::BasicCDRClusterSet():
	CacheableData()
{
	clusters_.clear();
	clusters_.resize(6, NULL);
}

BasicCDRClusterSet::BasicCDRClusterSet(utility::vector1<CDRClusterOP> const clusters):
	CacheableData()
{
	clusters_ = clusters;
}

BasicCDRClusterSet::BasicCDRClusterSet(const BasicCDRClusterSet& src):
	CacheableData()
{
	clusters_ = src.clusters_;
}

BasicCDRClusterSet::~BasicCDRClusterSet(){}

CacheableDataOP
BasicCDRClusterSet::clone() const{
	return CacheableDataOP( new BasicCDRClusterSet(*this) );
}

void
BasicCDRClusterSet::set_cluster( CDRNameEnum cdr, CDRClusterCOP cluster ){
	clusters_[ cdr ] = cluster->clone();
}	

void
BasicCDRClusterSet::set_clusters( utility::vector1<CDRClusterOP> const clusters ){
	assert( clusters.size() == 6 );
	clusters_ = clusters;
}

CDRClusterCOP
BasicCDRClusterSet::get_cluster(CDRNameEnum cdr) const {
	return clusters_[ cdr ];
}

//BasicCDRClusterSet::get_clusters() const {
//	return clusters_;
//}


} //clusters
} //antibody
} //protocols
