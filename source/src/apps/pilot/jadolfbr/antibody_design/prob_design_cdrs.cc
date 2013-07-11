// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jadolfbr/antibody_design/prob_design_cdrs.cc
/// @brief Designs CDRs using cluster probabilities.  Test App.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/design/AntibodyCDRDesigner.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

using namespace protocols::antibody;
using namespace protocols::antibody::design;

///@brief Wrapper to AntibodyCDRDesigner due to needed ab_info object.  Designs CDRs based on cluster using probability and conservative strategies.
///
class Designer : public protocols::moves::Mover {
	
	Designer(){};
	
	virtual ~Designer(){};
	
	virtual std::string
	get_name() const{
		return "AntibodyCDRDesignerWrapper";
	}
	
	void
	apply(core::pose::Pose pose){
		if (! protocols::antibody::check_if_pose_renumbered_for_clusters(pose)){
			utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.")
		}
		AntibodyInfoOP ab_info = new AntibodyInfo(pose, Modified_AHO);
		ab_info->show(std::cout);
		ab_info->setup_CDR_clusters(pose);
		AntibodyCDRDesigner des = new AntibodyCDRDesigner(ab_info);
		des->apply(pose);
	}
	
};




int main(int argc, char* argv[]){
	
	devel::init(argc, argv);

	
	try{
		protocols::jd2::JobDistributor::get_instance()->go(new Designer);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cout << "Exception: " << std::endl;
		excn.show( std::cerr );
	}
	
	return(0);
}

