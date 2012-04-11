// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/LHRepulsiveRamp.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/LHRepulsiveRamp.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>


#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>



#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::fmt;

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/antibody2/AntibodyUtil.hh>







#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.LHRepulsiveRamp");


using namespace core;
namespace protocols {
namespace antibody2 {




    
    
    
    
// default constructor
LHRepulsiveRamp::LHRepulsiveRamp() : Mover() {

}

    
LHRepulsiveRamp::LHRepulsiveRamp(loops::Loops loops_in ) : Mover() {
    user_defined_ = true;
    init(loops_in, false);
}
    
    
LHRepulsiveRamp::LHRepulsiveRamp(antibody2::AntibodyInfoOP antibody_in) : Mover() {
    user_defined_ = true;
    init(antibody_in->all_cdr_loops_,false);
}
    
LHRepulsiveRamp::LHRepulsiveRamp(antibody2::AntibodyInfoOP antibody_in, bool camelid) : Mover() {
    user_defined_ = true;
    init(antibody_in->all_cdr_loops_, camelid);
}
    
    
// default destructor
LHRepulsiveRamp::~LHRepulsiveRamp() {}
    
//clone
protocols::moves::MoverOP LHRepulsiveRamp::clone() const {
    return( new LHRepulsiveRamp() );
}
    
    

    
    
void LHRepulsiveRamp::init(loops::Loops loops_in, bool camelid ) 
{
    set_default();
    is_camelid_ = camelid;
    all_loops_ = loops_in;
    
    tf_ = new pack::task::TaskFactory;
    
    // score functions
    dock_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
    dock_scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
    dock_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. ); 
    
    pack_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
}
    
    
void LHRepulsiveRamp::set_default(){
    benchmark_= false;
    rep_ramp_cycles_= 3 ;
    min_type_="dfpmin_armijo_nonmonotone";
    rot_mag_ = 2.0 ;
    trans_mag_ = 0.1 ;
    temperature_ = 0.8;
    min_threshold_=15.0;
    cycles_= 4;
}
    
    
std::string LHRepulsiveRamp::get_name() const {
    return "LHRepulsiveRamp";
}

    
    
    
void LHRepulsiveRamp::finalize_setup(pose::Pose & pose ){
    setup_packer_task(pose, tf_);
    
    Size nres_ = pose.total_residue();
    
    ( *dock_scorefxn_ )( pose );
    
    //setting MoveMap
    cdr_dock_map_ = new kinematics::MoveMap();
    cdr_dock_map_->clear();
    cdr_dock_map_->set_chi( false );
    cdr_dock_map_->set_bb( false );
    utility::vector1< bool> bb_is_flexible( nres_, false );
    utility::vector1< bool> sc_is_flexible( nres_, false );

    select_loop_residues( pose, all_loops_, false/*include_neighbors*/, bb_is_flexible);
    cdr_dock_map_->set_bb( bb_is_flexible );
    select_loop_residues( pose, all_loops_, true/*include_neighbors*/, sc_is_flexible);
    cdr_dock_map_->set_chi( sc_is_flexible );
    cdr_dock_map_->set_jump( 1, true );
    for( Size ii = 2; ii <= all_loops_.num_loop() + 1; ii++ )
        cdr_dock_map_->set_jump( ii, false );

    
    
    //set up sidechain movers for rigid body jump and loop & neighbors
    utility::vector1_size rb_jump;
    rb_jump.push_back( 1 );
    using namespace core::pack::task;
    using namespace core::pack::task::operation;
    // selecting movable c-terminal residues
    ObjexxFCL::FArray1D_bool loop_residues( nres_, false );
    for( Size i = 1; i <= nres_; i++ ) {
        loop_residues(i) = sc_is_flexible[i]; 
    } // check mapping
    
    using namespace protocols::toolbox::task_operations;
    tf_->push_back( new RestrictToInterface( rb_jump, loop_residues ) );
}

    
    
    
///////////////////////////////////////////////////////////////////////////
/// @begin repulsive_ramp
///
/// @brief ramping up the fullatom repulsive weight slowly to allow the
///        partners to relieve clashes and make way for each other
///
/// @detailed This routine is specially targetted to the coupled
///           optimization of docking partners and the loop region.  The
///           loop modelling & all previous  steps  involve mainly
///           centroid  mode .On switching  on fullatom mode, one is bound
///           to end up with clashes.To relieve the clashes, it is
///           essential to slowly  dial up the  repulsive weight instead of
///           turning it on to the maximum  value in one single step
///
/// @param[in] input pose which is assumed to have a docking fold tree
///
/// @global_read fa_rep : fullatom repulsive weight
///
/// @global_write fa_rep ( It is reset to the original value at the end )
///
/// @remarks A particular portion is  commented out,which can be
///          uncommented if one  uses a  low  resolution  homology  model.
///          Check details in the beginning of the commented out region
///
/// @references
///
/// @authors Aroop 07/13/2010    
///
/// @last_modified 07/13/2010
///////////////////////////////////////////////////////////////////////////
    
void LHRepulsiveRamp::apply( pose::Pose & pose ) {
    
    
    finalize_setup(pose );


    
    
    
    // remove cutpoints variants for all cdrs
    // "true" forces removal of variants even from non-cutpoints
    loops::remove_cutpoint_variants( pose, true );
    
    using namespace core::chemical;
    for ( loops::Loops::const_iterator it = all_loops_.begin(),
         it_end = all_loops_.end();	it != it_end; ++it ) {
        core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, it->cut() );
        core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER,it->cut()+1);
    }
    // add scores to map
    ( *dock_scorefxn_ )( pose );
    
    // dampen fa_rep weight
    Real rep_weight_max = dock_scorefxn_->get_weight( core::scoring::fa_rep );
    Real func_tol(1.0);
    //mjo commenting out 'nb_list' because it is unused and causes a warning
    //bool nb_list( true );
    if( benchmark_ ) {
        rep_ramp_cycles_ = 1;
        cycles_ = 1;
        min_threshold_ = 150.0;
        func_tol = 10.0;
    }
    
    Real rep_ramp_step = (rep_weight_max - 0.02) / Real(rep_ramp_cycles_-1);
    for ( Size i = 1; i <= rep_ramp_cycles_; i++ ) {
        Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
        dock_scorefxn_->set_weight( core::scoring::fa_rep, rep_weight );
        snugfit_MC_min ( pose);
    }

    
    
}
    
    
    


    
    
    
    
    
    
void LHRepulsiveRamp::snugfit_MC_min ( pose::Pose & pose_in) {
        
		using namespace moves;
		bool nb_list = true;
        
        simple_moves::MinMoverOP min_mover = new simple_moves::MinMover( cdr_dock_map_, dock_scorefxn_,
                                                                        min_type_, min_threshold_, nb_list );
        
		//set up rigid body movers
        rigid::RigidBodyPerturbMoverOP rb_perturb=new rigid::RigidBodyPerturbMover(pose_in,
                                                                                   *cdr_dock_map_, rot_mag_, trans_mag_ , rigid::partner_downstream, true );
        

        
        simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( pack_scorefxn_, tf_ );
		SequenceMoverOP rb_mover = new SequenceMover;
		rb_mover->add_mover( rb_perturb );
		rb_mover->add_mover( pack_rottrial );
		JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_mover, min_mover, dock_scorefxn_,	min_threshold_ );
        

		MonteCarloOP mc = new MonteCarlo( pose_in, *dock_scorefxn_, temperature_ );
		TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc);
    
		RepeatMoverOP first_mcm_cycles = new RepeatMover( rb_mover_min_trial, cycles_ );
		first_mcm_cycles->apply( pose_in );
        
		return;
        
	} // snugfit_MC_min
    





    
    

    
    
    
    



} // namespace antibody2
} // namespace protocols





