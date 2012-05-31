// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineCDRH3HighRes.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/RefineCDRH3HighRes.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/H3RefineCCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>


#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below





using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.RefineCDRH3HighRes");




using namespace core;
namespace protocols {
namespace antibody2 {

    
    
    
// default constructor
RefineCDRH3HighRes::RefineCDRH3HighRes( ) : Mover() {}
    
RefineCDRH3HighRes::RefineCDRH3HighRes( AntibodyInfoOP antibody_info ) : Mover() {   
    user_defined_ = false;
    ab_info_ = antibody_info;
    init();
}
    

RefineCDRH3HighRes::RefineCDRH3HighRes( AntibodyInfoOP antibody_info, std::string refine_mode ) : Mover() {
    user_defined_ = false;
    ab_info_    = antibody_info;
    refine_mode_  = refine_mode;
    init();
}
  
    
RefineCDRH3HighRes::RefineCDRH3HighRes(AntibodyInfoOP antibody_info, 
                                       std::string refine_mode,
                                       scoring::ScoreFunctionCOP highres_scorefxn 
                                       ) : Mover()  {
    user_defined_ = true;
    ab_info_    = antibody_info;
    refine_mode_  = refine_mode;
    highres_scorefxn_ = new scoring::ScoreFunction(*highres_scorefxn);
    
    init();   
}
    

    
    
void RefineCDRH3HighRes::init( ) {
    set_default();
}

    
    
    
    
void RefineCDRH3HighRes::set_default()
{ 
    refine_mode_    = "legacy_refine_ccd";
    H3_filter_      = true;
    flank_relax_    = true;
    high_cst_       = 100.0;
    
    if(!user_defined_){
        highres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function("standard", "score12" );
            highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
            highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
            highres_scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );
    }
    
}
    
    
    
// default destructor
RefineCDRH3HighRes::~RefineCDRH3HighRes() {}
    
//clone
protocols::moves::MoverOP RefineCDRH3HighRes::clone() const {
    return( new RefineCDRH3HighRes() );
}
    
    

std::string RefineCDRH3HighRes::get_name() const {
    return "RefineCDRH3HighRes";
}

    
    
void RefineCDRH3HighRes::pass_start_pose(core::pose::Pose & start_pose){
    start_pose_ = start_pose;
}
    
    

    
    
void RefineCDRH3HighRes::apply(core::pose::Pose &pose){
    
    
    using namespace basic::options;
	using namespace basic::options::OptionKeys;
    
    
    if( refine_mode_ == "legacy_refine_ccd" ){
        H3RefineCCDOP legacy_refine_ccd = new H3RefineCCD(ab_info_, "h3", highres_scorefxn_); 
            legacy_refine_ccd -> pass_start_pose(start_pose_);
        legacy_refine_ccd -> apply(pose);
    }
    else{
        
        loops::LoopsOP pass_loops = new loops::Loops(); 
        pass_loops->add_loop(  *(ab_info_->get_CDR_loop("h3"))  );   
        
        if(refine_mode_ == "refine_ccd"){ 
            loops::loop_mover::refine::LoopMover_Refine_CCD refine_ccd( pass_loops, highres_scorefxn_ );
            if(get_native_pose()) {     refine_ccd.set_native_pose( get_native_pose() );    }
            refine_ccd.apply( pose );
        }
                             
        else if(refine_mode_ == "refine_kic"){
            //loops.remove_terminal_loops( pose );
            loops::loop_mover::refine::LoopMover_Refine_KIC refine_kic( pass_loops, highres_scorefxn_ );
            if(get_native_pose()) {     refine_kic.set_native_pose( get_native_pose() );    }
            refine_kic.apply( pose );

        }
        else{
            utility_exit_with_message("the refinement method is not available!");
        }
    }
    
    
}
    
    





} // namespace antibody2
} // namespace protocols



