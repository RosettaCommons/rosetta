// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.cc
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


#include <utility/exit.hh>


static basic::Tracer TR("protocols.antibody.design.AntibodyDesignMover");
namespace protocols{
namespace antibody {
namespace design {
	
using namespace basic::options;
using namespace basic::options::OptionKeys::antibody;
using namespace core::scoring;
using namespace protocols::antibody;
using namespace core::import_pose;
using core::pose::Pose;
using std::string;

AntibodyDesignMover::AntibodyDesignMover() :
		protocols::moves::Mover(),
		//ab_info_(ab_info),
		//cdr_db_parser_(new AntibodyDatabaseManager()),
		//ab_info_(NULL),
		graft_designer_(NULL),
		modeler_(NULL),
		scorefxn_(NULL)
{
	ab_info_ = NULL;
	read_options();
	cdr_db_parser_ = new AntibodyDatabaseManager();
	protocols::moves::Mover::type( "AntibodyDesign" );
}

    
AntibodyDesignMover::~AntibodyDesignMover(){}
    
std::string
AntibodyDesignMover::get_name() const {
	return "AntibodyDesignMover";
}
    
protocols::moves::MoverOP
AntibodyDesignMover::clone() const {
	return new AntibodyDesignMover(*this);
}

void
AntibodyDesignMover::read_options(){
	return;
}

void
AntibodyDesignMover::set_graft_designer(AntibodyGraftDesignerOP & graft_designer){
	graft_designer_=graft_designer;
}

void
AntibodyDesignMover::run_graft_designer(core::pose::Pose& pose){
	if (! graft_designer_){
		graft_designer_ = new AntibodyGraftDesigner(ab_info_);
	}
	
	graft_designer_->apply(pose);
}

void
AntibodyDesignMover::model_starting_pose(core::pose::Pose& pose){
    	//Probably better off relaxing with coordinate constraints..
	//bool relax_cdrs = option[design::modeler::relax_starting_cdrs]();
	//if (relax_cdrs){
	//TR<<"Minimizing CDRs"<<std::endl;
	protocols::antibody::set_harmonic_constraints(ab_info_, pose);
	//modeler_->relax_cdrs(pose);
        }

void
AntibodyDesignMover::set_ab_info(AntibodyInfoOP & ab_info){
	ab_info_ = ab_info;
}

void
AntibodyDesignMover::set_scorefxn(ScoreFunctionOP & scorefxn){
	scorefxn_ = scorefxn;
}

void
AntibodyDesignMover::apply(core::pose::Pose& pose){
	
	///Setup Objects///
	if (! ab_info_){
		ab_info_ = new AntibodyInfo(pose, Modified_AHO);
		ab_info_->setup_CDR_clusters(pose);
		ab_info_->show(TR);

	}
	
	if (ab_info_->get_Current_AntibodyNumberingScheme()!="Modifed_AHO"){
		utility_exit_with_message("Antibody Design Protocol requires the Modified AHO numbering scheme");
	}
	
	TR << "Mover not yet complete - Lets check if pose is loaded." <<std::endl;
	if (! scorefxn_){
		scorefxn_ = getScoreFunction();
	}
        
	(*scorefxn_)(pose);
	scorefxn_->show(TR, pose);
        
	if (!modeler_){
		modeler_ = new AntibodyDesignModeler(ab_info_);
	}
	
	TR << "Initialize Options" << std::endl;

	//model_starting_pose(pose);
	
	TR <<"Running Graft Design" <<std::endl;
	
	run_graft_designer(pose);
	
	TR <<"Exiting..."<< std::endl;

}
    
    
} //Design
} //Antibody
} //Protocols
