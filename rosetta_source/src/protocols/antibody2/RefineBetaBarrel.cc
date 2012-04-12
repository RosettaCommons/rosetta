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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/tools/make_vector1.hh>



#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.antibody2.RefineBetaBarrel");

using namespace core;

namespace protocols {
namespace antibody2 {
        
        
        
        
// default constructor
RefineBetaBarrel::RefineBetaBarrel() : Mover() {}
RefineBetaBarrel::~RefineBetaBarrel() {}
    

RefineBetaBarrel::RefineBetaBarrel(AntibodyInfoOP antibody_info) : Mover() {        
    user_defined_ = false;
    ab_info_ = antibody_info;

    init();
}
    
RefineBetaBarrel::RefineBetaBarrel(AntibodyInfoOP antibody_info,
                                    core::scoring::ScoreFunctionCOP dock_scorefxn,
                                    core::scoring::ScoreFunctionCOP pack_scorefxn) : Mover() {
    user_defined_ = true;
    ab_info_ = antibody_info;
    dock_scorefxn_ = new core::scoring::ScoreFunction(*dock_scorefxn);
    pack_scorefxn_ = new core::scoring::ScoreFunction(*pack_scorefxn);
    
    init();
}
    

void RefineBetaBarrel::init( ){
    
    if(!user_defined_){
        dock_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
            dock_scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
            dock_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
        pack_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
    }

    repulsive_ramp_ = true;

    // set up objects
    lh_repulsive_ramp_ = new LHRepulsiveRamp(ab_info_, dock_scorefxn_, pack_scorefxn_);
    dock_mcm_protocol_ = new docking::DockMCMProtocol( ab_info_->LH_dock_jump(), dock_scorefxn_, pack_scorefxn_ );
        //TODO: check move map, fold tree, task factory that you should pass into.
        //lh_snugfit_ = new LHSnugFitLegacy(ab_info_);

}
   

    

    
    
std::string RefineBetaBarrel::get_name() const {
    return "RefineBetaBarrel";
}
    
    
    
    
    
    
void RefineBetaBarrel::apply( pose::Pose & pose ) {

    all_cdr_VL_VH_fold_tree( pose, ab_info_->all_cdr_loops_ );
    

    if(repulsive_ramp_) {
        lh_repulsive_ramp_->apply(pose);
    }
    
    TR<<"finish repulsive ramping !"<<std::endl;

    
    //lh_snugfit_ ->set_task_factory(tf_);
    //lh_snugfit_ -> apply(pose);
    
    // TODO: 
    // JQX: check fold tree, move map, and task factory
    dock_mcm_protocol_ -> apply(pose);

    TR<<"finish dock_mcm ramping !"<<std::endl;
}


    
    


    
    
    

}// namespace antibody2
}// namespace protocols





