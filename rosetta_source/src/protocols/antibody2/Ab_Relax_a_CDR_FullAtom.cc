// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_Relax_a_CDR_FullAtom.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/Ab_Relax_a_CDR_FullAtom.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


#include <protocols/antibody2/AntibodyUtil.hh>

#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below





using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.Ab_Relax_a_CDR_FullAtom");




using namespace core;
namespace protocols {
namespace antibody2 {

    
    
    
// default constructor
Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( ) : Mover() 
{
    set_default();
    
    is_camelid_ = false;
    init();
}

Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( AntibodyInfoOP antibody_info, std::string loop_name ) : Mover() 
{
    set_default();

    is_camelid_ = false;
    ab_info_    = antibody_info;
    loop_name_  = loop_name;
    
    
    init();
}
  
    

    
Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( bool is_camelid, AntibodyInfoOP antibody_info ) : Mover() 
{        
    set_default();
        

    is_camelid_ = is_camelid;
    ab_info_ = antibody_info;
    init();
}
    
    

    
void Ab_Relax_a_CDR_FullAtom::set_default(){ 
    benchmark_ = false;
    include_neighbors_ = true;
    max_cycle_close_trial_ = 20;
	refine_input_loop_ = false;
    flank_relax_ = true;
	flank_size_ = 2;
	h3_random_cut_ = false;
	H3_filter_ = true;
    min_tolerance_ = 0.001;
    high_move_temp_ = 2.00;
    minimization_type_ = "dfpmin_armijo_nonmonotone" ;
    init_temp_ = 2.0;
    last_temp_ = 0.5;
    gamma_ = std::pow( (last_temp_/init_temp_), (1.0/inner_cycles_));
    use_pymol_diy_ = false;
    neighbor_dist_ = 10.0;
}
    
    
    
// default destructor
Ab_Relax_a_CDR_FullAtom::~Ab_Relax_a_CDR_FullAtom() {}
    
//clone
protocols::moves::MoverOP Ab_Relax_a_CDR_FullAtom::clone() const {
    return( new Ab_Relax_a_CDR_FullAtom() );
}
    
    

    
    
void Ab_Relax_a_CDR_FullAtom::init( ) 
{

    the_loop_   = *(ab_info_->get_CDR_loop(loop_name_));
    loop_begin_ = the_loop_.start();
    loop_end_   = the_loop_.stop();
    loop_size_  = ( loop_end_ - loop_begin_ ) + 1;
    
    cutpoint_ = loop_begin_ + int(loop_size_/2);

    if( h3_random_cut_ ){
        //	cutpoint_ = dle_choose_random_cutpoint(loop_begin_, loop_end_);
    }

    the_loop_.set_cut(cutpoint_);
    
    n_small_moves_ =  numeric::max(Size(5), Size(loop_size_/2)) ;
    inner_cycles_ = loop_size_;
    outer_cycles_ = 2; //JQX: assume the SnugFit step has done some minimization
    if(  refine_input_loop_ ){
        outer_cycles_ = 5;
    }
    if( benchmark_ ) {
        min_tolerance_ = 1.0;
        n_small_moves_ = 1;
        inner_cycles_ = 1;
        outer_cycles_ = 1;
    }
    
    highres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function("standard", "score12" );
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
	// adding constraints
	//highres_scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );


    
}
    

    
    
std::string Ab_Relax_a_CDR_FullAtom::get_name() const {
    return "Ab_Relax_a_CDR_FullAtom";
}

    
    
void Ab_Relax_a_CDR_FullAtom::pass_start_pose(core::pose::Pose & start_pose){
    start_pose_ = start_pose;
}
    
    
    
    
    
    
    
    
    
void Ab_Relax_a_CDR_FullAtom::finalize_setup( core::pose::Pose & pose ){

    setup_packer_task( start_pose_, tf_ );
    
    // set cutpoint variants for correct chainbreak scoring
    if( !pose.residue( cutpoint_ ).is_upper_terminus() ) {
        if( !pose.residue( cutpoint_ ).has_variant_type(chemical::CUTPOINT_LOWER)){
            core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint_ );
        }
        if( !pose.residue( cutpoint_ + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) ){
            core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint_ + 1 );
        }
    }

    // the list of residues that are allowed to pack
    for(Size ii=1; ii <=pose.total_residue();ii++) {allow_repack_.push_back(false);}
    select_loop_residues( pose, the_loop_, include_neighbors_, allow_repack_, neighbor_dist_);

    
    // the list of residues that are allowed to change backbone
    utility::vector1< bool> allow_bb_move( pose.total_residue(), false );
    for( Size ii = loop_begin_; ii <= loop_end_; ii++ ){
        allow_bb_move[ ii ] = true;
    }
    
    // the movemap of the h3 loop, if flank_relax_=false, flank_cdrh3_map_=cdrh3_map_
    cdrh3_map_ = new kinematics::MoveMap();
    cdrh3_map_->set_jump( 1, false );
    cdrh3_map_->set_bb( allow_bb_move );
    cdrh3_map_->set_chi( allow_repack_ );
    
    if( flank_relax_) {
        utility::vector1< bool> flank_allow_bb_move( allow_bb_move );
        for( Size i = 1; i <= pose.total_residue(); i++ ){
            if(  (i >= (loop_begin_ - flank_size_)) && (i <= (loop_end_ + flank_size_))   ){
                flank_allow_bb_move[i] = true;
            }
        }
       
        flank_cdrh3_map_ = new kinematics::MoveMap();
        flank_cdrh3_map_->set_jump( 1, false );
        flank_cdrh3_map_->set_bb( flank_allow_bb_move );
        flank_cdrh3_map_->set_chi( allow_repack_ );
    }
    else{
        flank_cdrh3_map_ = cdrh3_map_;
    }
    

    // below are the definitions of a bunch of movers 
    using namespace protocols;
    using namespace protocols::simple_moves;
    using namespace protocols::loops;
    using namespace protocols::moves;
    using namespace protocols::toolbox::task_operations;
    using namespace pack;
    using namespace pack::task;
    using namespace pack::task::operation;
    using loop_closure::ccd::CcdMover;
    using loop_closure::ccd::CcdMoverOP;
    
    
    // the Monte Carlo mover
    mc_ = new protocols::moves::MonteCarlo( pose, *highres_scorefxn_, init_temp_ );
    
    
    // setup some easy objects to change the fold trees
    // if not flank_relax =false, then   change_FT_to_flankloop_ = change_FT_to_simpleloop_
    simple_fold_tree( pose, loop_begin_ - 1, cutpoint_, loop_end_ + 1 );
    change_FT_to_simpleloop_ = new ChangeFoldTreeMover( pose.fold_tree() );
    
    if(flank_relax_){
        simple_fold_tree( pose, loop_begin_ - flank_size_ - 1, cutpoint_, loop_end_ + flank_size_ + 1 );
    }
    change_FT_to_flankloop_ = new ChangeFoldTreeMover( pose.fold_tree() );
    
    
    // pack the loop and its neighboring residues
    loop_repack_ = new PackRotamersMover(highres_scorefxn_);
    setup_packer_task( start_pose_, tf_ );
    ( *highres_scorefxn_ )( pose );
    tf_->push_back( new RestrictToInterface( allow_repack_ ) );
    loop_repack_->task_factory(tf_);
        //loop_repack_->apply( pose_in );
    
    
    // minimize amplitude of moves if correct parameter is set
    BackboneMoverOP small_mover = new SmallMover( cdrh3_map_, high_move_temp_, n_small_moves_ );
    BackboneMoverOP shear_mover = new ShearMover( cdrh3_map_, high_move_temp_, n_small_moves_ );

    small_mover->angle_max( 'H', 2.0 ); small_mover->angle_max( 'E', 5.0 ); small_mover->angle_max( 'L', 6.0 );
    shear_mover->angle_max( 'H', 2.0 ); shear_mover->angle_max( 'E', 5.0 ); shear_mover->angle_max( 'L', 6.0 );

    
    // ccd moves
    CcdMoverOP ccd_moves = new CcdMover( the_loop_, cdrh3_map_ );
    RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves_);
    
    
    // minimization mover
    loop_min_mover_ = new MinMover( flank_cdrh3_map_, highres_scorefxn_, minimization_type_, min_tolerance_, true /*nb_list*/ );

    
    // put everything into a sequence mover
    wiggle_cdr_h3_ = new SequenceMover() ;
    wiggle_cdr_h3_->add_mover( change_FT_to_simpleloop_ );
    wiggle_cdr_h3_->add_mover( small_mover );     if(use_pymol_diy_) wiggle_cdr_h3_->add_mover(pymol_);
    wiggle_cdr_h3_->add_mover( shear_mover );     if(use_pymol_diy_) wiggle_cdr_h3_->add_mover(pymol_);
    wiggle_cdr_h3_->add_mover( ccd_cycle );       if(use_pymol_diy_) wiggle_cdr_h3_->add_mover(pymol_);
    wiggle_cdr_h3_->add_mover( change_FT_to_flankloop_ );
    wiggle_cdr_h3_->add_mover( loop_min_mover_ ); if(use_pymol_diy_) wiggle_cdr_h3_->add_mover(pymol_);

}
    

    
    
    
    
    
    
    
    
    
    
    
//APPLY
void Ab_Relax_a_CDR_FullAtom::apply( pose::Pose & pose ) {
    using namespace protocols::simple_moves;
    using namespace protocols::moves;
    using namespace protocols::toolbox::task_operations;
    using namespace pack::task;
    using namespace pack::task::operation;
        
    
    if ( !pose.is_fullatom() ){utility_exit_with_message("Fullatom poses only");}
    
    TR <<  " Relaxing Fullatom CDR H3 loop" << std::endl;
    

    finalize_setup(pose);

    
    bool closed_cutpoints( false );
    Size cycle( 1 );
    while( !closed_cutpoints && cycle<max_cycle_close_trial_ ) {
            TR <<  "    Refining CDR H3 loop in HighRes.  close_trial_cycle="<<cycle << std::endl;

        
            change_FT_to_flankloop_->apply(pose);
            loop_min_mover_->apply(pose);
        
            // rotamer trials
            select_loop_residues( pose, the_loop_, include_neighbors_, allow_repack_, neighbor_dist_);
            setup_packer_task( start_pose_, tf_ );
            ( *highres_scorefxn_ )( pose );
            tf_->push_back( new RestrictToInterface( allow_repack_ ) );
            RotamerTrialsMoverOP pack_rottrial = new RotamerTrialsMover( highres_scorefxn_, tf_ );
            pack_rottrial->apply( pose );        
        
        
        

            Real temperature = init_temp_;
            mc_->reset( pose ); // monte carlo reset
        
            bool relaxed_H3_found_ever( false );
            if( H3_filter_){
                relaxed_H3_found_ever = CDR_H3_filter( pose,the_loop_, is_camelid_);
            }
        
            // outer cycle
            for(Size i = 1; i <= outer_cycles_; i++) {
                mc_->recover_low( pose );
                Size h3_attempts(0); // number of H3 checks after refinement
            
                // inner cycle
                for ( Size j = 1; j <= inner_cycles_; j++ ) {
                    mc_->set_temperature( temperature*= gamma_ );
                
                    wiggle_cdr_h3_->apply( pose );
                
                    // rotamer trials
                    select_loop_residues( pose, the_loop_, include_neighbors_, allow_repack_, neighbor_dist_);
                    setup_packer_task( start_pose_, tf_ );
                    ( *highres_scorefxn_ )( pose );
                    tf_->push_back( new RestrictToInterface( allow_repack_ ) );
                    pack_rottrial->task_factory(tf_);
                    pack_rottrial->apply( pose );
                
                    bool relaxed_H3_found_current(false);
                    // H3 filter check
                    if(H3_filter_ && (h3_attempts <= inner_cycles_)) {
                        h3_attempts++;
                        relaxed_H3_found_current = CDR_H3_filter(pose, the_loop_, is_camelid_);
                    
                        if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
                            mc_->boltzmann( pose );
                        }
                        else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
                            relaxed_H3_found_ever = true;
                            mc_->reset( pose );
                        }
                        else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
                            --j;
                            continue;
                        }
                        else if( relaxed_H3_found_ever && relaxed_H3_found_current ){
                            mc_->boltzmann( pose );
                        }
                    }
                    else {
                        if( H3_filter_ ) {
                            bool relaxed_H3_found_current(false);
                            relaxed_H3_found_current = CDR_H3_filter(pose,the_loop_, is_camelid_);
                            if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
                                mc_->boltzmann( pose );
                            }
                            else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
                                relaxed_H3_found_ever = true;
                                mc_->reset( pose );
                            }
                            else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
                                mc_->recover_low( pose );
                            }
                            else if( relaxed_H3_found_ever && relaxed_H3_found_current ){
                                mc_->boltzmann( pose );
                            }
                        }
                        else{
                            mc_->boltzmann( pose );
                        }
                    }
                
                    if ( numeric::mod(j,Size(20))==0 || j==inner_cycles_ ) {
                        // repack trial
                        loop_repack_ = new PackRotamersMover( highres_scorefxn_ );
                        setup_packer_task( start_pose_, tf_ );
                        ( *highres_scorefxn_ )( pose );
                        tf_->push_back( new RestrictToInterface( allow_repack_ ) );
                        loop_repack_->task_factory( tf_ );
                        loop_repack_->apply( pose );
                        mc_->boltzmann( pose );
                    }
                
                
                } // inner 
            } // outer
            mc_->recover_low( pose );

            // minimize
            change_FT_to_flankloop_->apply( pose );
            loop_min_mover_->apply( pose );
        



            closed_cutpoints = cutpoints_separation( pose, ab_info_ );
            ++cycle;
    } // while( ( cut_separation > 1.9 )
        
    TR << "Finished Relaxing CDR H3 Loop" << std::endl;
        
    return;
    
} 
    

    

    




} // namespace antibody2
} // namespace protocols






















/*
 
 JQX: Just record what I deleted
 
 
 
 
 void CDRH3Modeler2::build_fullatom_loop( core::pose::Pose & pose ) {
    using namespace core::pose;
    using namespace core::scoring;
    using namespace protocols::moves;
 
    if( !apply_fullatom_mode_ )
        return;
 
    TR <<  "H3M Modeling Fullatom CDR H3 loop" << std::endl;
 
    antibody2::AntibodyInfo starting_antibody;
    starting_antibody = ab_info_;
    bool closed_cutpoints( false );
 
    Size cycle( 1 );
    while( !closed_cutpoints && cycle<max_cycle_close_trial_ ) {
        ab_info_ = starting_antibody;
        loop_fa_relax( pose, ab_info_.get_CDR_loop("h3")->start(),
        ab_info_.get_CDR_loop("h3")->stop()-1 + base_ );
        closed_cutpoints = cutpoints_separation( pose, ab_info_ );
        ++cycle;
    } // while( ( cut_separation > 1.9 )
 
    TR <<  "H3M Finished modeling Fullatom CDR H3 loop" << std::endl;
 
    return;
 } // build_fullatom_loop
 
 
 
 
 if( is_camelid_ ) {
    bool store_current_loop = current_loop_is_H3_;
    bool store_H3_filter = H3_filter_;
    current_loop_is_H3_ = false;
    H3_filter_ = false;
 
    antibody2::AntibodyInfo starting_antibody;
    starting_antibody = ab_info_;
    bool closed_cutpoints( false );
 
    Size cycle (1 );
    while( !closed_cutpoints && cycle < max_cycle_close_trial_ ) {
        ab_info_ = starting_antibody;
        loop_fa_relax( pose_in, ab_info_.get_CDR_loop("h2")->start(),
        ab_info_.get_CDR_loop("h2")->stop()  );
        closed_cutpoints = cutpoints_separation( pose_in, ab_info_ );
        ++cycle;
    } // while( ( cut_separation > 1.9 )
 
    // Restoring variables to initial state
    current_loop_is_H3_ = store_current_loop;
    H3_filter_ = store_H3_filter;
}

 
 */
