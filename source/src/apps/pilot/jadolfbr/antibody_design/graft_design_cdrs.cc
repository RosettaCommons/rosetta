// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jadolfbr/antibody_design/graft_design_cdrs.cc
/// @brief Designs an antibody using GraftDesign. 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>

#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <devel/init.hh>




using namespace protocols::antibody;
using namespace protocols::antibody::design;

//Documentation:  This app designs antibodies using North cluster graft design. 
// Not recommended for general use.  Though it is being tested for designs in the wet lab. 
// If you would like to use this app, please email me for help.  It requires a separate AntibodyDatabase, which is not yet downloadable.

class GraftDesignCDRs : public protocols::moves::Mover{
public:
	GraftDesignCDRs(){};
	
	virtual ~GraftDesignCDRs(){};
	
	virtual
	std::string
	get_name() const {
		return "GraftDesignCDRs";
	}
	
	//virtual
	//protocols::moves::MoverOP
	//fresh_instance() const {
	//	return new GraftDesignCDRs;;
	//}
	
	void
	apply(core::pose::Pose & pose){
		AntibodyInfoOP ab_info = new AntibodyInfo(pose, Modified_AHO);
		ab_info->show(std::cout);
		ab_info->setup_CDR_clusters(pose);
	
		AntibodyGraftDesignerOP ab_designer = new AntibodyGraftDesigner(ab_info);
		ab_designer->apply(pose);
		utility::vector1< core::pose::PoseOP > result_poses;
		result_poses = ab_designer->get_top_designs();
		
		protocols::jd2::JobOP current_job( protocols::jd2::JobDistributor::get_instance()->current_job());
		
		for (core::Size i=2; i<=result_poses.size(); ++i){
			std::string tag = "graft_ensemble_"+utility::to_string(i)+"_";
			//Need to have it obey path_.  Make this change AFTER tests.
			
			std::cout << "Outputting ensemble " << i << std::endl;
			
			//Add information for other poses to be output.
			//ab_info->setup_CDR_clusters(*(result_poses[i]));
			//for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
			//	CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			//	std::pair<CDRClusterEnum, core::Real> result = ab_info->get_CDR_cluster(cdr_name);
			//	std::string output = "REMARK CLUSTER "+ ab_info->get_cluster_name(result.first) +" "+utility::to_string(result.second);
			//std::cout << output;
			//	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(output);
			//}
			
			protocols::jd2::JobDistributor::get_instance()->job_outputter()->other_pose(current_job, *(result_poses[i]), tag);
		}
		
		//Add information to final pose
		//ab_info->setup_CDR_clusters(pose);
		for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
			CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			std::pair<CDRClusterEnum, core::Real> result = ab_info->get_CDR_cluster(cdr_name);
			std::string output = "REMARK CLUSTER "+ ab_info->get_cluster_name(result.first) +" "+utility::to_string(result.second);
			//std::cout << output <<std::endl;
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(output);
		}
	}
};

int main(int argc, char* argv[]){
	
	try{

		devel::init(argc, argv);

	
		protocols::jd2::JobDistributor::get_instance()->go(new GraftDesignCDRs);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cout << "Exception: " << std::endl;
		excn.show( std::cerr );
	}
	
	return(0);
}
