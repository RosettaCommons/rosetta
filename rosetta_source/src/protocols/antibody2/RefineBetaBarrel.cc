// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineBetaBarrel.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/RefineBetaBarrel.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/LHRepulsiveRamp.hh>
#include <protocols/antibody2/LHSnugFitLegacy.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/docking/DockMCMProtocol.hh>

#include <utility/tools/make_vector1.hh>



#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.antibody2.RefineBetaBarrel");

using namespace core;

namespace protocols {
namespace antibody2 {
        
        
        
        
// default constructor
RefineBetaBarrel::RefineBetaBarrel() : Mover() {
    user_defined_ = false;
}

RefineBetaBarrel::RefineBetaBarrel(AntibodyInfoOP antibody_info) : Mover() {        
    user_defined_ = true;
    init(antibody_info);
}
    

void RefineBetaBarrel::init(AntibodyInfoOP antibody_info){
    repulsive_ramp_ = true;
    ab_info_ = antibody_info;
    lh_repulsive_ramp_ = new LHRepulsiveRamp(ab_info_);
    lh_snugfit_ = new LHSnugFitLegacy(ab_info_);
    //TODO: replace the vector1 by ab_info->LH_jump;
    dock_mcm_protocol_ = new docking::DockMCMProtocol();

}
   
    
RefineBetaBarrel::~RefineBetaBarrel() {}
    
    
void RefineBetaBarrel::set_default(){
}
    
    
std::string RefineBetaBarrel::get_name() const {
    return "RefineBetaBarrel";
}
    
    
    
void RefineBetaBarrel::apply( pose::Pose & pose ) {

    all_cdr_VL_VH_fold_tree( pose, ab_info_->all_cdr_loops_ );
    

    if(repulsive_ramp_) {lh_repulsive_ramp_->apply(pose);}
    
    TR<<"finish repulsive ramping"<<std::endl;

    
    //lh_snugfit_ ->set_task_factory(tf_);
    //lh_snugfit_ -> apply(pose);
    
    // TODO: 
    // JQX: check fold tree, move map, and task factory
    dock_mcm_protocol_ -> apply(pose);
    
}


    
    


    
    
    

}// namespace antibody2
}// namespace protocols





