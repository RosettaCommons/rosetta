// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyDesignModeler.cc
/// @brief Handles modeling of the antibody.  Before and after design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Antibody Headers
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

// Core Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

//Protocol Headers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/simple_moves/MinMover.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>

static basic::Tracer TR("protocols.antibody.design.AntibodyDesignModeler");
namespace protocols {
namespace antibody {
namespace design {
            
            
using namespace core::kinematics;         
using namespace core::pose;    
using namespace core::scoring;      
using namespace basic::options;      
using namespace basic::options::OptionKeys::antibody;        
using namespace protocols::simple_moves;                      
using core::Size;      
using std::string;
            
        
AntibodyDesignModeler::AntibodyDesignModeler(AntibodyInfoOP & ab_info){
	ab_info_ = ab_info;
	cdrs_.resize(CDRNameEnum_total);
	set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
}
        
        
AntibodyDesignModeler::~AntibodyDesignModeler(){}
        

void
AntibodyDesignModeler::set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting){
	for (core::Size i=cdr_start; i<=cdr_end; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		set_cdr(cdr, setting);
	}
}

void
AntibodyDesignModeler::set_cdr_only(CDRNameEnum cdr, bool setting){
	if (setting==true){
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, false);
		set_cdr(cdr, true);
	}
	else{
		set_cdr_range(CDRNameEnum_start, CDRNameEnum_total, true);
		set_cdr(cdr, false);
	}
}
void
AntibodyDesignModeler::set_cdr(CDRNameEnum const cdr, bool setting){
	cdrs_[cdr]=setting;
}

//void
//AntibodyDesignModeler::copy_cdrs(core::pose::Pose & from_pose, core::pose::Pose & to_pose, core::Size overhang){
	
//}

void
AntibodyDesignModeler::relax_cdrs(Pose & pose, ScoreFunctionOP const & scorefxn, bool centroid_mode){

	TR << "Beginning to refine structure..." << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn->clone();
	
	MoveMapOP mm = new MoveMap();
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		if (cdrs_[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			//Add any options here!
			ab_info_->add_CDR_to_MoveMap(pose, mm, cdr);
		}
	}

	if (centroid_mode){
		protocols::relax::CentroidRelaxOP rel = new protocols::relax::CentroidRelax(mm);
		rel->apply(pose);
		return;
	}
	else{
		local_scorefxn->set_weight(chainbreak, 100.00); //First PyRosetta Toolkit lesson.
		scorefxn->show(pose);
		protocols::relax::FastRelaxOP rel = new protocols::relax::FastRelax(scorefxn);//Stack allocation and destroy.
		rel->set_movemap(mm);
		rel->apply(pose);
		scorefxn->show(pose);
	}
}

void
AntibodyDesignModeler::minimize_cdrs(Pose & pose, ScoreFunctionOP const & scorefxn){
	TR << "Beginning to refine structure..." << std::endl;
	ScoreFunctionOP local_scorefxn = scorefxn->clone();
	MoveMapOP mm = new MoveMap();
	for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		ab_info_->add_CDR_to_MoveMap(pose, mm, cdr);
	}
	scorefxn->show(pose);
	protocols::simple_moves::MinMover min_mover(mm, local_scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false /*use_nblist*/ );
	min_mover.apply(pose);
	scorefxn->show(pose);
}

//void
//AntibodyDesignModeler::repack_CDRs(Pose & pose, ScoreFunctionOP const scorefxn){
	
//}
 
//void
//AntibodyDesignModeler::run_snugdock(Pose & pose, ScoreFunctionOP const scorefxn){
	
//}


}//design
}//antibody
}//protocols
