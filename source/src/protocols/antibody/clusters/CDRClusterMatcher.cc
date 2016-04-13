// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterMatcher.cc
/// @brief Simple class for identifying CDR cluster in an antibody
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/clusters/CDRClusterMatcher.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/clusters/util.hh>

#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/string_util.hh>

#include <numeric/NumericTraits.hh>


#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/types.hh>
#include <math.h>
#include <map>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

// Boost headers
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.cluster.CDRClusterMatcher" );


namespace protocols {
namespace antibody {
namespace clusters {
using namespace protocols::antibody;
using utility::vector1;
using core::Real;

CDRClusterMatcher::CDRClusterMatcher(){
	using namespace basic::options;

	center_cluster_db_path_="sampling/antibodies/cluster_center_dihedrals.txt";
	
	allow_rama_mismatches_ = option [ OptionKeys::antibody::allow_omega_mismatches_for_north_clusters]();
	load_center_data();
	
}

CDRClusterMatcher::~CDRClusterMatcher(){}

void
CDRClusterMatcher::load_center_data(){

	//I need to turn this text file into a database soon.
	CDRClusterEnumManagerOP cluster_man( new CDRClusterEnumManager() );
	AntibodyEnumManagerCOP ab_man( AntibodyEnumManagerOP( new AntibodyEnumManager() ) );



	std::string cdr_name;
	std::string cluster;
	std::string fullcluster_name;
	std::string phis;
	std::string psis;
	std::string omegas;
	std::string cis_trans_conf;
	core::Size length;
	std::string type;

	utility::io::izstream clus_stream;
	basic::database::open(clus_stream, center_cluster_db_path_);
	Size line_count = 0;
	while ( ! clus_stream.eof() ) {
		++line_count;
		clus_stream >> cdr_name >> length >> cluster >> type >> fullcluster_name >> cis_trans_conf >> phis >> psis >> omegas;
		struct ClusterData single_cluster_data;
		vector1< std::string > phiSP;
		vector1< std::string > psiSP;
		single_cluster_data.cdr = ab_man->cdr_name_string_to_enum(cdr_name);
		single_cluster_data.cluster = cluster_man->cdr_cluster_string_to_enum(fullcluster_name);
		single_cluster_data.length = length;
		single_cluster_data.cis_trans_conf = cis_trans_conf;

		boost::split(phiSP, phis, boost::is_any_of(","));
		boost::split(psiSP, psis, boost::is_any_of(","));

		for ( core::Size i = 1; i <= phiSP.size(); ++i ) {
			std::istringstream phi_stream(phiSP[i]);
			std::istringstream psi_stream(psiSP[i]);
			Real phi; phi_stream >> phi;
			Real psi; psi_stream >> psi;
			single_cluster_data.phis.push_back(phi);
			single_cluster_data.psis.push_back(psi);
		}
		cluster_data_.push_back(single_cluster_data);
	}
	clus_stream.close();
}

CDRClusterOP
CDRClusterMatcher::get_cdr_cluster(core::pose::Pose const & pose, CDRNameEnum const cdr, core::Size const start, core::Size const end) const {

	std::string cis_trans_conf= get_pose_cis_trans_conformation(pose, start, end);

	std::map< std::string, vector1< Real > > pose_angles = get_pose_angles(pose, start, end);
	std::map <Real, ClusterData > k_distances_to_cluster;
	vector1< Real > k_distances;

	core::Size length = end - start +1;
	bool cluster_found = false;
	TR << "Length: " << length << " Omega: "<< cis_trans_conf << std::endl;
	bool cis_trans_match = false;
	bool length_match = false;

	for ( core::Size i = 1; i <= cluster_data_.size(); ++i ) {
		ClusterData data = cluster_data_[i];
		
		//Overall Booleans
		if ( (data.cdr ==  cdr) && (data.length ==  length) ) length_match = true;
		if (  (data.cdr ==  cdr) && (data.length ==  length) && data.cis_trans_conf == cis_trans_conf ) cis_trans_match = true;
		
		
		if ( data.cdr ==  cdr && data.length ==  length && data.cis_trans_conf == cis_trans_conf) {
			cluster_found = true;

			core::Real k_distance_to_cluster = calculate_dihedral_distance(data.phis, pose_angles["phi"], data.psis, pose_angles["psi"]);
			k_distances_to_cluster[k_distance_to_cluster] = data;
			k_distances.push_back(k_distance_to_cluster);
		} else { continue; }
	}
	
	if ( ( ! cluster_found ) && length_match && allow_rama_mismatches_){
		for ( core::Size i = 1; i <= cluster_data_.size(); ++i ) {
			ClusterData data = cluster_data_[i];

			if ( (data.cdr ==  cdr) && (data.length ==  length) ) {
				cluster_found = true;
				core::Real k_distance_to_cluster = calculate_dihedral_distance(data.phis, pose_angles["phi"], data.psis, pose_angles["psi"]);
				k_distances_to_cluster[k_distance_to_cluster] = data;
				k_distances.push_back(k_distance_to_cluster);
			} else { continue; }
		}
	}

	///Take the minimum distance as the cluster.
	CDRClusterEnum cluster;
	core::Real distance;
	if ( ! cluster_found ) {
		cluster = NA;
		distance = 1000;
		
		if (! length_match ){
			TR.Warning << std::endl << TR.Red << "No known cluster of CDR length "<< length << "found. Setting as NA" << TR.Reset << std::endl;
		}
		else if ( ! cis_trans_match ){
			TR.Warning << std::endl << TR.Red << "*** No known cluster of CDR length " << length << " omega of " << cis_trans_conf << " found. ***" << TR.Reset << std::endl;
			TR.Warning << TR.Red << "*** Consider using the command-line option -allow_omega_mismatches_for_north_clusters to find the closest cluster! ***" << TR.Reset << std::endl << std::endl;
		}
		
	} else {

		//Get minimum and set cluster.
		distance = utility::min(k_distances);
		cluster = k_distances_to_cluster[distance].cluster;

	}

	CDRClusterOP cdr_cluster( new CDRCluster(pose, cdr, length, cluster, start, distance, cis_trans_match) );
	return cdr_cluster;


}

CDRClusterOP
CDRClusterMatcher::get_closest_cluster(core::pose::Pose const & pose, core::Size const start, core::Size const end) const {
	std::string cis_trans_conf= get_pose_cis_trans_conformation(pose, start, end);

	std::map< std::string, vector1< Real > > pose_angles = get_pose_angles(pose, start, end);
	std::map <Real, ClusterData> k_distances_to_cluster;
	vector1< Real > k_distances;

	core::Size length = end - start +1;
	bool cluster_found = false;
	bool length_match = false;
	bool cis_trans_match = false;
	
	
	for ( core::Size i = 1; i <= cluster_data_.size(); ++i ) {
		ClusterData data = cluster_data_[i];
		
		if ( data.length == length) length_match = true;
		if ( data.length == length && data.cis_trans_conf == cis_trans_conf ) cis_trans_match = true;
		
		if ( data.length == length && data.cis_trans_conf == cis_trans_conf ) {
			cluster_found = true;
			
			core::Real k_distance_to_cluster = calculate_dihedral_distance(data.phis, pose_angles["phi"], data.psis, pose_angles["psi"]);
			k_distances_to_cluster[k_distance_to_cluster] = data;
			k_distances.push_back(k_distance_to_cluster);
		} else { continue; }
	}

	if ( (! cluster_found ) && length_match && allow_rama_mismatches_){
		for ( core::Size i = 1; i <= cluster_data_.size(); ++i ) {
			ClusterData data = cluster_data_[i];

			if ( data.length ==  length) {
				cluster_found = true;
				core::Real k_distance_to_cluster = calculate_dihedral_distance(data.phis, pose_angles["phi"], data.psis, pose_angles["psi"]);
				k_distances_to_cluster[k_distance_to_cluster] = data;
				k_distances.push_back(k_distance_to_cluster);
			} else { continue; }
		}
	}
	
	
	///Take the minimum distance as the cluster.

	CDRClusterEnum cluster;
	core::Real distance;
	CDRNameEnum cdr;
	if ( ! cluster_found ) {
		cluster = NA;
		distance = 1000;
		cdr = l1;
		if (! length_match ){
			TR.Warning << std::endl << TR.Red << "No known cluster of length "<< length << "found. Setting as NA" << TR.Reset << std::endl;
		}
		else if ( ! cis_trans_match ){
			TR.Warning << std::endl << TR.Red << "*** No known cluster of length " << length << " and omega of " << cis_trans_conf << " found. ***" << TR.Reset << std::endl;
			TR.Warning << TR.Red << "*** Consider using the command-line option -allow_omega_mismatches_for_north_clusters to find the closest cluster! ***" << TR.Reset << std::endl << std::endl;
		}
	} else {

		//Get minimum and set cluster.
		distance = utility::min(k_distances);
		cluster = k_distances_to_cluster[distance].cluster;
		cdr = k_distances_to_cluster[distance].cdr;

	}

	CDRClusterOP cdr_cluster( new CDRCluster(pose, cdr, length, cluster, start, distance) );
	return cdr_cluster;

}

std::map< std::string, vector1< core::Real > >
CDRClusterMatcher::get_pose_angles(core::pose::Pose const & pose, core::Size const start, core::Size const end) const {

	std::map< std::string, vector1< core::Real > > pose_angles;
	for ( Size resnum = start; resnum<=end; resnum++ ) {
		pose_angles["phi"].push_back(pose.phi(resnum));
		pose_angles["psi"].push_back(pose.psi(resnum));
	}
	return pose_angles;
}


} //clusters
} //antibody
}//protocols
