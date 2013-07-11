// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jadolfbr/cluster_utilities/identify_cdr_clusters.cc
///
/// @brief This identifies CDR clusters of antibody.  Modified_AHO numbering required.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <devel/init.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

//Options
//#include <basic/options/option.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/string_util.hh>

using namespace protocols::antibody;


//Documentation:  This application identifies the CDR cluster in an antibody, renumbered with Modified_AHO (Used by North clusters).  Works with one PDB.  Prints to screen
//  Use: Renumber antibody using [http://dunbrack.fccc.edu/IgClassify/] (Not quite done).  Outputs info, and appends it to a new PDB that it will write.  
//  Reference: North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256.
class IdentifyCDRClusters : public protocols::moves::Mover{
public:
	IdentifyCDRClusters(){};
	
	virtual ~IdentifyCDRClusters(){};
	
	virtual
	std::string
	get_name() const {
		return "IdentifyCDRClusters";
	}
	
	void
	apply(core::pose::Pose & pose){
		
		if (! protocols::antibody::check_if_pose_renumbered_for_clusters(pose)){
			utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please visit www.xxx.edu");
		}
		AntibodyInfoOP ab_info = new AntibodyInfo(pose, Modified_AHO);
		ab_info->show(std::cout);
		ab_info->setup_CDR_clusters(pose);
		
		
		for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
			CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			std::pair<CDRClusterEnum, core::Real> result = ab_info->get_CDR_cluster(cdr_name);
			std::string output = "REMARK CLUSTER "+ ab_info->get_cluster_name(result.first) +" "+utility::to_string(result.second);
			//std::cout << output;
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(output);
		}
		//std::string reference = "REF: North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256.";
		//protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(reference);
	}
};

int main(int argc, char* argv[])
{
	devel::init(argc, argv);
	
	try{
		protocols::jd2::JobDistributor::get_instance()->go(new IdentifyCDRClusters);
		std::cout<< "Please reference North, B., A. Lehmann, and R.L. Dunbrack, Jr., A new clustering of antibody CDR loop conformations. JMB, 2011. 406(2): p. 228-56."<<std::endl;
	} catch(utility::excn::EXCN_Base & excn){
		std::cout << "Exception"<<std::endl;
		excn.show(std::cerr);
	}
	
	return(0);
}

//To be included:
// 1) Multiple PDB's (done through JD2)
// 2) Output to an optional text file (done through additional info at end of PDB)
// 3) Optionally Output to a file instead of PDB
