// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/clusters/util.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/util.hh>

#include <numeric/NumericTraits.hh>
#include <core/pose/util.hh>
#include <core/pose/PDB_Info.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <math.h>
#include <cmath>


namespace protocols {
namespace antibody {
namespace clusters {
	using utility::vector1;
	using core::Real;
	using core::Size;

static basic::Tracer TR("protocols.antibody.clusters");

std::string
get_pose_cis_trans_conformation(core::pose::Pose const & pose, core::Size const start, core::Size const end){
	std::string model_ss;
	core::Real cis_cutoff = 90.0;
	for (core::Size resnum = start-1; resnum <= end-1; ++resnum) {
		if (std::abs(pose.omega(resnum)) >= cis_cutoff) {
			model_ss = model_ss+"T";
		} else {
			model_ss = model_ss+"C";
		}
	}
	return model_ss;
}

core::Real
calculate_dihedral_distance(vector1< Real > cluster_phis, vector1< Real> pose_phis, vector1< Real > cluster_psis, vector1< Real > pose_psis){
	Real k_distance_to_cluster = 0.0;
	Real PI = numeric::NumericTraits<Real>::pi();
	for (Size i=1; i<=cluster_phis.size(); ++i) {

		Real phi_d = (2 * (1- cos ((pose_phis[i]-cluster_phis[i])*PI/180)));
		Real psi_d = (2 * (1- cos ((pose_psis[i]-cluster_psis[i])*PI/180)));
		k_distance_to_cluster = k_distance_to_cluster+phi_d+psi_d;

	}
	return k_distance_to_cluster;
}

bool
check_if_pose_renumbered_for_clusters(core::pose::Pose const & pose){
	

	if (core::pose::has_chain("L", pose)){
		//L1
		if (pose.residue(pose.pdb_info()->pdb2pose('L', 23)).name1() != 'C'){
			TR << "Problem with L1 starting anchor residue name" <<std::endl;
			return false;
		}
		if (pose.residue(pose.pdb_info()->pdb2pose('L', 43)).name1() != 'W'){
			TR << "Problem with L1 ending anchor residue name" <<std::endl;
			return false;
		}
		
		//L3
		if (pose.residue(pose.pdb_info()->pdb2pose('L', 106)).name1() != 'C'){
			TR << "Problem with L3 starting anchor residue name" <<std::endl;
			return false;
		}
		if (pose.residue(pose.pdb_info()->pdb2pose('L', 139)).name1() != 'F'){
			TR << "Problem with L3 ending anchor residue name" <<std::endl;
			return false;
		}
	}
	
	if (core::pose::has_chain("H", pose)){
		//H1
		if (pose.residue(pose.pdb_info()->pdb2pose('H', 23)).name1() != 'C'){
			TR << "Problem with H1 starting anchor residue name" <<std::endl;
			return false;
		}
		if (pose.residue(pose.pdb_info()->pdb2pose('H', 43)).name1() != 'W'){
			TR << "Problem with H1 ending anchor residue name" <<std::endl;
			return false;
		}
		
		//H3
		if (pose.residue(pose.pdb_info()->pdb2pose('H', 106)).name1() != 'C'){
			TR << "Problem with H3 starting anchor residue name" <<std::endl;
			return false;
		}
		if (pose.residue(pose.pdb_info()->pdb2pose('H', 139)).name1() != 'W'){
			TR << "Problem with H3 ending anchor residue name" <<std::endl;
			return false;
		}
	}
	
	return true;
}

} //clusters
} //antibody
} //protocols
