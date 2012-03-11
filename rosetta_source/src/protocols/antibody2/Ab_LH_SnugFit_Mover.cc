// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_LH_SnugFit_Mover.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/Ab_LH_SnugFit_Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody2/Ab_Info.hh>
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







#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.Ab_LH_SnugFit_Mover");
using namespace core;





namespace protocols {
namespace antibody2 {




    
    
    
    
// default constructor
Ab_LH_SnugFit_Mover::Ab_LH_SnugFit_Mover() : Mover() {

}

    
Ab_LH_SnugFit_Mover::Ab_LH_SnugFit_Mover(loops::Loops loops_in ) : Mover() {
    user_defined_ = true;
    init(loops_in, false);
}
    
    
Ab_LH_SnugFit_Mover::Ab_LH_SnugFit_Mover(antibody2::Ab_Info & antibody_in) : Mover() {
    user_defined_ = true;
    init(antibody_in.all_cdr_loops_,false);
}
    
Ab_LH_SnugFit_Mover::Ab_LH_SnugFit_Mover(antibody2::Ab_Info & antibody_in, bool camelid) : Mover() {
    user_defined_ = true;
    init(antibody_in.all_cdr_loops_, camelid);
}
    
    
// default destructor
Ab_LH_SnugFit_Mover::~Ab_LH_SnugFit_Mover() {}
    
//clone
protocols::moves::MoverOP Ab_LH_SnugFit_Mover::clone() const {
    return( new Ab_LH_SnugFit_Mover() );
}
    
    

    
    
void Ab_LH_SnugFit_Mover::init(loops::Loops loops_in, bool camelid ) 
{
    is_camelid_ = camelid;
    all_loops_ = loops_in;
}
    
    
void Ab_LH_SnugFit_Mover::set_default(){
        
}
    
    
std::string Ab_LH_SnugFit_Mover::get_name() const {
    return "Ab_LH_SnugFit_Mover";
}

    
    
    
    
    
    
    
    
    
//APPLY
void Ab_LH_SnugFit_Mover::apply( pose::Pose & pose ) {

    
    snugfit_mcm_protocol( pose, all_loops_ ) ;
    
    
}
    
    
    




void Ab_LH_SnugFit_Mover::snugfit_mcm_protocol( pose::Pose & pose_in, loops::Loops loops_in ) 
{
    
    using namespace moves;
    bool nb_list = true;
    Size nres = pose_in.total_residue();
    
    //MC move
    Real trans_mag ( 0.1 );
    Real rot_mag ( 5.0 );
    
    // rb minimization
    std::string min_type = "dfpmin_armijo_nonmonotone";
    Real min_threshold ( 15.0 ); /* score unit */
    
    // score functions
    using namespace core::scoring;
    core::scoring::ScoreFunctionOP scorefxn;
    scorefxn = core::scoring::ScoreFunctionFactory::
        create_score_function( "docking", "docking_min" );
    scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
    scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );
    
    // score functions
    core::scoring::ScoreFunctionOP pack_scorefxn;
    pack_scorefxn = core::scoring::ScoreFunctionFactory::
    create_score_function( "standard" );
    
    // remove cutpoints variants for all cdrs
    // "true" forces removal of variants even from non-cutpoints
    loops::remove_cutpoint_variants( pose_in, true );
    
    using namespace core::chemical;
    for ( loops::Loops::const_iterator it = loops_in.begin(),
         it_end = loops_in.end();	it != it_end; ++it ) {
        core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_LOWER, it->cut() );
        core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_UPPER,it->cut()+1);
    }
    
    //setting MoveMap
    kinematics::MoveMapOP cdr_dock_map;
    cdr_dock_map = new kinematics::MoveMap();
    cdr_dock_map->clear();
    cdr_dock_map->set_chi( false );
    cdr_dock_map->set_bb( false );
    utility::vector1< bool> is_flexible( nres, false );
    bool include_neighbors( false );
    select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
    cdr_dock_map->set_bb( is_flexible );
    include_neighbors = true;
    select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
    cdr_dock_map->set_chi( is_flexible );
    cdr_dock_map->set_jump( 1, true );
    for( Size ii = 2; ii <= loops_in.num_loop() + 1; ii++ )
        cdr_dock_map->set_jump( ii, false );
    
    
    //set up minimizer movers
    simple_moves::MinMoverOP min_mover = new simple_moves::MinMover( cdr_dock_map, scorefxn, min_type, min_threshold, nb_list );
    
    //set up rigid body movers
    rigid::RigidBodyPerturbMoverOP rb_perturb = new rigid::RigidBodyPerturbMover( pose_in,
                                                                                 *cdr_dock_map, rot_mag, trans_mag, rigid::partner_downstream, true );
    
    setup_packer_task( pose_in );
    
    //set up sidechain movers for rigid body jump and loop & neighbors
    utility::vector1_size rb_jump;
    rb_jump.push_back( 1 );
    using namespace core::pack::task;
    using namespace core::pack::task::operation;
    
    // selecting movable c-terminal residues
    ObjexxFCL::FArray1D_bool loop_residues( nres, false );
    for( Size i = 1; i <= nres; i++ )
        loop_residues( i ) = is_flexible[ i ]; // check mapping
    using namespace protocols::toolbox::task_operations;
    tf_->push_back( new RestrictToInterface( rb_jump, loop_residues ) );
    
    
    
    simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( pack_scorefxn, tf_ );
    
    simple_moves::PackRotamersMoverOP pack_interface_repack = new simple_moves::PackRotamersMover( pack_scorefxn );
    pack_interface_repack->task_factory(tf_);
    
    Real temperature = 0.8;
    MonteCarloOP mc = new MonteCarlo( pose_in, *scorefxn, temperature );
    
    TrialMoverOP pack_interface_trial = new TrialMover(pack_interface_repack, mc );
    
    protocols::docking::SidechainMinMoverOP scmin_mover = new
    protocols::docking::SidechainMinMover( core::scoring::ScoreFunctionCOP( pack_scorefxn ), core::pack::task::TaskFactoryCOP( tf_ ) );
    TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc );
    
    SequenceMoverOP rb_mover = new SequenceMover;
    rb_mover->add_mover( rb_perturb );
    rb_mover->add_mover( pack_rottrial );
    
    JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_mover, min_mover, scorefxn, min_threshold);
    TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc  );
    
    SequenceMoverOP repack_step = new SequenceMover;
    repack_step->add_mover( rb_mover_min_trial );
    repack_step->add_mover( pack_interface_trial );
    repack_step->add_mover( scmin_trial );
    
    CycleMoverOP rb_mover_min_trial_repack  = new CycleMover;
    for ( Size i=1; i < 8; ++i )
        rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
    rb_mover_min_trial_repack->add_mover( repack_step );
    
    //set up initial repack mover
    SequenceMoverOP initial_repack = new SequenceMover;
    initial_repack->add_mover( pack_interface_trial );
    initial_repack->add_mover( scmin_trial );
    
    //set up initial and final min_trial movers for docking
    TrialMoverOP minimize_trial = new TrialMover( min_mover, mc );
    
    //set up mcm cycles and mcm_repack cycles
    RepeatMoverOP mcm_four_cycles = new RepeatMover( rb_mover_min_trial, 4 );
    
    Size cycles = 3;
    if ( benchmark_ ) cycles = 1;
    RepeatMoverOP mcm_final_cycles = new RepeatMover( rb_mover_min_trial_repack, cycles );
    
    SequenceMoverOP snugfit_mcm = new SequenceMover;
    snugfit_mcm->add_mover( initial_repack );
    snugfit_mcm->add_mover( minimize_trial );
    snugfit_mcm->add_mover( mcm_four_cycles );
    snugfit_mcm->add_mover( mcm_final_cycles );
    snugfit_mcm->add_mover( minimize_trial );
    
    snugfit_mcm->apply ( pose_in );
    
    return;
} // snugfit_mcm_protocol


    
    
    
    //TODO:
    //JQX:
    // you saw the functino of setup_packer_task was used multiple times
    // in relaxCDR and repulsive_ramp as well, it is better to make a separate utility 
    // class, or put it in the ab_info_ object
void Ab_LH_SnugFit_Mover::setup_packer_task(pose::Pose & pose_in ) {
    
    using namespace pack::task;
    using namespace pack::task::operation;
        
    if( init_task_factory_ ) {
        tf_ = new TaskFactory( *init_task_factory_ );
        TR << "AbModeler Reinitializing Packer Task" << std::endl;
        return;
    }
    else
        tf_ = new TaskFactory;
        
    TR << "AbModeler Setting Up Packer Task" << std::endl;
        
    tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
    tf_->push_back( new InitializeFromCommandline );
    tf_->push_back( new IncludeCurrent );
    tf_->push_back( new RestrictToRepacking );
    tf_->push_back( new NoRepackDisulfides );
        
    // incorporating Ian's UnboundRotamer operation.
    // note that nothing happens if unboundrot option is inactive!
    pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new pack::rotamer_set::UnboundRotamersOperation();
    unboundrot->initialize_from_command_line();
    operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
    tf_->push_back( unboundrot_operation );
    // adds scoring bonuses for the "unbound" rotamers, if any
    core::pack::dunbrack::load_unboundrot( pose_in );
        
    init_task_factory_ = tf_;
        
    TR << "AbModeler Done: Setting Up Packer Task" << std::endl;
        
} // setup_packer_task
    
    
    
    
    



} // namespace antibody2
} // namespace protocols





