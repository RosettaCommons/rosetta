// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_LH_RepulsiveRamp_Mover.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/Ab_LH_RepulsiveRamp_Mover.hh>

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
#include <protocols/antibody2/Ab_util.hh>







#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below



using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.Ab_LH_RepulsiveRamp_Mover");
using namespace core;





namespace protocols {
namespace antibody2 {




    
    
    
    
// default constructor
Ab_LH_RepulsiveRamp_Mover::Ab_LH_RepulsiveRamp_Mover() : Mover() {

}

    
Ab_LH_RepulsiveRamp_Mover::Ab_LH_RepulsiveRamp_Mover(loops::Loops loops_in ) : Mover() {
    user_defined_ = true;
    init(loops_in, false);
}
    
    
Ab_LH_RepulsiveRamp_Mover::Ab_LH_RepulsiveRamp_Mover(antibody2::Ab_Info & antibody_in) : Mover() {
    user_defined_ = true;
    init(antibody_in.all_cdr_loops_,false);
}
    
Ab_LH_RepulsiveRamp_Mover::Ab_LH_RepulsiveRamp_Mover(antibody2::Ab_Info & antibody_in, bool camelid) : Mover() {
    user_defined_ = true;
    init(antibody_in.all_cdr_loops_, camelid);
}
    
    
// default destructor
Ab_LH_RepulsiveRamp_Mover::~Ab_LH_RepulsiveRamp_Mover() {}
    
//clone
protocols::moves::MoverOP Ab_LH_RepulsiveRamp_Mover::clone() const {
    return( new Ab_LH_RepulsiveRamp_Mover() );
}
    
    

    
    
void Ab_LH_RepulsiveRamp_Mover::init(loops::Loops loops_in, bool camelid ) 
{
    is_camelid_ = camelid;
    all_loops_ = loops_in;
}
    
    
void Ab_LH_RepulsiveRamp_Mover::set_default(){
        
}
    
    
std::string Ab_LH_RepulsiveRamp_Mover::get_name() const {
    return "Ab_LH_RepulsiveRamp_Mover";
}

    
    
    
    
    
    
    
    
    
//APPLY
void Ab_LH_RepulsiveRamp_Mover::apply( pose::Pose & pose ) {

    
    repulsive_ramp( pose, all_loops_ ) ;
    
    
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
void Ab_LH_RepulsiveRamp_Mover::repulsive_ramp(
                                  pose::Pose & pose_in,
                                  loops::Loops loops_in ) {
        
		Size nres = pose_in.total_residue();
        
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
        
        
		// score functions
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
		// add scores to map
		( *scorefxn )( pose_in );
        
		// dampen fa_rep weight
		Real rep_weight_max = scorefxn->get_weight( core::scoring::fa_rep );
		Size rep_ramp_cycles(3);
		Size cycles(4);
		Real minimization_threshold(15.0);
		Real func_tol(1.0);
		//mjo commenting out 'nb_list' because it is unused and causes a warning
		//bool nb_list( true );
		if( benchmark_ ) {
			rep_ramp_cycles = 1;
			cycles = 1;
			minimization_threshold = 150.0;
			func_tol = 10.0;
		}
        
		Real rep_ramp_step = (rep_weight_max - 0.02) / Real(rep_ramp_cycles-1);
		for ( Size i = 1; i <= rep_ramp_cycles; i++ ) {
			Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
			scorefxn->set_weight( core::scoring::fa_rep, rep_weight );
			snugfit_MC_min ( pose_in, cdr_dock_map, cycles, minimization_threshold,
                            scorefxn, pack_scorefxn, is_flexible);
		}
        
		return;
	} // repulsive_ramp
    
void Ab_LH_RepulsiveRamp_Mover::snugfit_MC_min (
                                   pose::Pose & pose_in,
                                   kinematics::MoveMapOP cdr_dock_map,
                                   Size cycles,
                                   Real minimization_threshold,
                                   core::scoring::ScoreFunctionOP scorefxn,
                                   core::scoring::ScoreFunctionOP pack_scorefxn,
                                   utility::vector1< bool> is_flexible ) {
        
		using namespace moves;
		bool nb_list = true;
		Size nres = pose_in.total_residue();
        
        simple_moves::MinMoverOP min_mover = new simple_moves::MinMover( cdr_dock_map, scorefxn,
                                                                        "dfpmin_armijo_nonmonotone", minimization_threshold, nb_list );
        
		//set up rigid body movers
        rigid::RigidBodyPerturbMoverOP rb_perturb=new rigid::RigidBodyPerturbMover(pose_in,
                                                                                   *cdr_dock_map, 2.0, 0.1 , rigid::partner_downstream, true );
        
		setup_packer_task( pose_in, tf_ );
		//set up sidechain movers for rigid body jump and loop & neighbors
		utility::vector1_size rb_jump;
		rb_jump.push_back( 1 );
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		// selecting movable c-terminal residues
		ObjexxFCL::FArray1D_bool loop_residues( nres, false );
		for( Size i = 1; i <= nres; i++ ) loop_residues( i ) = is_flexible[ i ]; // check mapping
		using namespace protocols::toolbox::task_operations;
		tf_->push_back( new RestrictToInterface( rb_jump, loop_residues ) );
        
        simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( pack_scorefxn, tf_ );
		SequenceMoverOP rb_mover = new SequenceMover;
		rb_mover->add_mover( rb_perturb );
		rb_mover->add_mover( pack_rottrial );
		JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_mover, min_mover, scorefxn,	minimization_threshold );
        
		Real temperature = 0.8;
		MonteCarloOP mc = new MonteCarlo( pose_in, *scorefxn, temperature );
		TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc);
		RepeatMoverOP first_mcm_cycles = new RepeatMover( rb_mover_min_trial, cycles );
		first_mcm_cycles->apply( pose_in );
        
		return;
        
	} // snugfit_MC_min
    





    
    

    
    
    
    



} // namespace antibody2
} // namespace protocols





