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

#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>

#include <protocols/simple_moves/FragmentMover.hh>
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
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


#include <protocols/antibody2/Ab_util.hh>

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
Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom() : Mover() {

}

   
    
Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( bool current_loop_is_H3, bool H3_filter ) : Mover() 
{
    set_default();
    
    current_loop_is_H3_ = current_loop_is_H3;
    H3_filter_ = H3_filter;
    is_camelid_ = false;
    init();
}

Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( bool current_loop_is_H3, bool H3_filter, Ab_InfoOP antibody_info ) : Mover() 
{
    set_default();
        
    current_loop_is_H3_ = current_loop_is_H3;
    H3_filter_ = H3_filter;
    is_camelid_ = false;
    ab_info_ = antibody_info;
    init();
}
  
    
Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( bool current_loop_is_H3, bool H3_filter, bool is_camelid ) : Mover() 
{        
    set_default();
    
    current_loop_is_H3_ = current_loop_is_H3;
    H3_filter_ = H3_filter;
    is_camelid_ = is_camelid;
    init();
}
    
    Ab_Relax_a_CDR_FullAtom::Ab_Relax_a_CDR_FullAtom( bool current_loop_is_H3, bool H3_filter, bool is_camelid, Ab_InfoOP antibody_info ) : Mover() 
    {        
        set_default();
        
        current_loop_is_H3_ = current_loop_is_H3;
        H3_filter_ = H3_filter;
        is_camelid_ = is_camelid;
        ab_info_ = antibody_info;
        init();
    }
    
void Ab_Relax_a_CDR_FullAtom::set_default(){ 
    apply_fullatom_mode_ = true;
    max_cycle_ = 20;
    base_ = 1;
	refine_input_loop_ = true;
	snug_fit_ = true;
    freeze_h3_ = true;
    flank_relax_ = true;
	h3_flank_ = 2;
	decoy_loop_cutpoint_ = 0;
	h3_random_cut_ = false;
	min_base_relax_ = false;
    antibody_refine_ = true;

	antibody_build_ = true;
	H3_filter_ = true;
    current_loop_is_H3_ = true;


}
    
    
    
// default destructor
Ab_Relax_a_CDR_FullAtom::~Ab_Relax_a_CDR_FullAtom() {}
    
//clone
protocols::moves::MoverOP Ab_Relax_a_CDR_FullAtom::clone() const {
    return( new Ab_Relax_a_CDR_FullAtom() );
}
    
    

    
    
void Ab_Relax_a_CDR_FullAtom::init( ) 
{
    if(is_camelid_) {
        snug_fit_=false;
    }

}
    
    
void Ab_Relax_a_CDR_FullAtom::setup_objects(){
        
}
    
    
std::string Ab_Relax_a_CDR_FullAtom::get_name() const {
    return "Ab_Relax_a_CDR_FullAtom";
}

    
    
void Ab_Relax_a_CDR_FullAtom::pass_start_pose(core::pose::Pose & start_pose){
    start_pose_ = start_pose;
}
    
    
    
    
    
    
//APPLY
void Ab_Relax_a_CDR_FullAtom::apply( pose::Pose & pose ) {


    build_fullatom_loop(pose)    ;
    
}
    
    

    void Ab_Relax_a_CDR_FullAtom::build_fullatom_loop( core::pose::Pose & pose ) {
        using namespace core::pose;
        using namespace core::scoring;
        using namespace protocols::moves;
        
        if( !apply_fullatom_mode_ )
            return;
        
        TR <<  "H3M Modeling Fullatom CDR H3 loop" << std::endl;
        
        antibody2::Ab_InfoOP starting_antibody;
        starting_antibody = ab_info_;
        bool closed_cutpoints( false );
        
        Size cycle( 1 );
        while( !closed_cutpoints && cycle<max_cycle_ ) {
            ab_info_ = starting_antibody;
            if(current_loop_is_H3_ && H3_filter_) {
                loop_fa_relax( pose, ab_info_->get_CDR_loop("h3")->start(), 
                                    ab_info_->get_CDR_loop("h3")->stop()-1 + base_ );
            }
            else if(!current_loop_is_H3_ && !H3_filter_&&is_camelid_){
                loop_fa_relax( pose, ab_info_->get_CDR_loop("h2")->start(),
                              ab_info_->get_CDR_loop("h2")->stop()  );
            }

            closed_cutpoints = cutpoints_separation( pose, ab_info_ );
            ++cycle;
        } // while( ( cut_separation > 1.9 )
        
        TR <<  "H3M Finished modeling Fullatom CDR H3 loop" << std::endl;
        
        return;
    } // build_fullatom_loop
    

    

    
    ///////////////////////////////////////////////////////////////////////////
    /// @begin loop_fa_relax  //JQX: this is actually the H3_refinment
    ///
    /// @brief actually relaxes the region specified
    ///
    /// @detailed This is all done in high resolution.Hence there are no rigid
    ///           body moves relative to the docking partners. Only small moves
    ///           are carried out here to see if there are better fits.
    ///           Repacking is carried out extensively after each move.
    ///
    /// @param[in] pose, loop begin position, loop end position
    ///
    /// @global_read none
    ///
    /// @global_write none
    ///
    /// @remarks
    ///
    /// @references
    ///
    /// @authors Aroop 02/04/2010
    ///
    /// @last_modified 02/04/2010
    ///////////////////////////////////////////////////////////////////////////
    void Ab_Relax_a_CDR_FullAtom::loop_fa_relax(
                                      pose::Pose & pose_in,
                                      Size const loop_begin,
                                      Size const loop_end )
    {
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
        
        TR << "H3M Relaxing CDR H3 Loop" << std::endl;
        
        // storing starting fold tree
        kinematics::FoldTree tree_in( pose_in.fold_tree() );
        
        //setting MoveMap
        kinematics::MoveMapOP cdrh3_map;
        cdrh3_map = new kinematics::MoveMap();
        cdrh3_map->clear();
        cdrh3_map->set_chi( false );
        cdrh3_map->set_bb( false );
        utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
        for( Size ii = loop_begin; ii <= loop_end; ii++ )
            allow_bb_move[ ii ] = true;
        cdrh3_map->set_bb( allow_bb_move );
        cdrh3_map->set_jump( 1, false );
        
        
        // minimize_set_local_min( false, 0 );// all non-move-list rsds minimized
        
        Size loop_size = ( loop_end - loop_begin ) + 1;
        Size cutpoint = loop_begin + int(loop_size/2);
        if( current_loop_is_H3_ ) {
            if( (antibody_build_ || antibody_refine_ ) &&
               !min_base_relax_ && !h3_random_cut_ &&
               (decoy_loop_cutpoint_ != 0))
                cutpoint = decoy_loop_cutpoint_;
            //else if( h3_random_cut_ )
            //	cutpoint = dle_choose_random_cutpoint(loop_begin, loop_end);
        }
        else
            cutpoint = loop_begin + Size( loop_size / 2);
        /*
         if( snug_fit_ && loops_flag_ && docking_local_refine_ &&
         dle_flag_ ) {
         one_loop = dle_ns::dle_loops;
         }
         else */
        loops::Loop one_loop( loop_begin, loop_end,	cutpoint,	0, false );
        
        // sets up stuff to use rosetta's fullatom energy function
        //initialize_fullatom();
        // maximum number of rotamers to allow (buried, surface)
        //set_rot_limit( 45, 27 );
        // maximum sum for the dunbrak prob (buried,surf)
        //set_perc_limit( 0.9999, 0.9999 );
        // include extra rotamers in chi1/chi2/chi1aro
        //design::active_rotamer_options.set_ex12( true, true, true);
        // checking if old rosetta full atom flag is on
        
        if ( !pose_in.is_fullatom() )
            utility_exit_with_message("Fullatom poses only");
        
        ChangeFoldTreeMoverOP one_loop_fold_tree;
        ChangeFoldTreeMoverOP with_flank_fold_tree;
        simple_fold_tree( pose_in, loop_begin - 1, cutpoint, loop_end + 1 );
        one_loop_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );
        with_flank_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );
        
        //////////////////
        // setup fold_tree
        utility::vector1< bool> flank_allow_bb_move( allow_bb_move  );
        if( current_loop_is_H3_  && flank_relax_ && freeze_h3_) {
            simple_fold_tree( pose_in, loop_begin - h3_flank_ - 1, cutpoint,
                             loop_end + h3_flank_ + 1 );
            with_flank_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );
            for( Size i = 1; i <= pose_in.total_residue(); i++ )
                if(  (i >= (loop_begin - h3_flank_)) && (i <= (loop_end + h3_flank_))   )
                    flank_allow_bb_move[i] = true;
        }
        else
            one_loop_fold_tree->apply( pose_in );
        
        // set cutpoint variants for correct chainbreak scoring
        if( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
            if( !pose_in.residue( cutpoint ).has_variant_type(chemical::CUTPOINT_LOWER))
                core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, cutpoint );
            if( !pose_in.residue( cutpoint + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
                core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, cutpoint + 1 );
        }
        
        
        
        utility::vector1< bool> allow_repack( pose_in.total_residue(), false );
        select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
                             allow_repack);
        cdrh3_map->set_chi( allow_repack );
        
        PackRotamersMoverOP loop_repack=new PackRotamersMover(highres_scorefxn_);
        setup_packer_task( start_pose_, tf_ );
        ( *highres_scorefxn_ )( pose_in );
        tf_->push_back( new RestrictToInterface( allow_repack ) );
        loop_repack->task_factory(tf_);
        // loop_repack->apply( pose_in );
        
        Real min_tolerance = 0.001;
        if( benchmark_ ) min_tolerance = 1.0;
        std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
        bool nb_list = true;
        MinMoverOP loop_min_mover = new MinMover( cdrh3_map,
                                                 highres_scorefxn_, min_type, min_tolerance, nb_list );
        
        // more params
        Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
        Size inner_cycles( loop_size );
        Size outer_cycles( 1 );
        if( antibody_refine_ || refine_input_loop_ )
            outer_cycles = 5;
        if( antibody_refine_ && snug_fit_ )
            outer_cycles = 2;
        if( benchmark_ ) {
            n_small_moves = 1;
            inner_cycles = 1;
            outer_cycles = 1;
        }
        
        Real high_move_temp = 2.00;
        // minimize amplitude of moves if correct parameter is set
        BackboneMoverOP small_mover = new SmallMover( cdrh3_map, high_move_temp, n_small_moves );
        BackboneMoverOP shear_mover = new ShearMover( cdrh3_map, high_move_temp, n_small_moves );
        if( min_base_relax_ ) {
            small_mover->angle_max( 'H', 0.5 );
            small_mover->angle_max( 'E', 0.5 );
            small_mover->angle_max( 'L', 1.0 );
            shear_mover->angle_max( 'H', 0.5 );
            shear_mover->angle_max( 'E', 0.5 );
            shear_mover->angle_max( 'L', 1.0 );
        }
        else {
            small_mover->angle_max( 'H', 2.0 );
            small_mover->angle_max( 'E', 5.0 );
            small_mover->angle_max( 'L', 6.0 );
            shear_mover->angle_max( 'H', 2.0 );
            shear_mover->angle_max( 'E', 5.0 );
            shear_mover->angle_max( 'L', 6.0 );
        }
        
        CcdMoverOP ccd_moves = new CcdMover( one_loop, cdrh3_map );
        RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);
        
        SequenceMoverOP wiggle_cdr_h3( new SequenceMover() );
        wiggle_cdr_h3->add_mover( one_loop_fold_tree );
        wiggle_cdr_h3->add_mover( small_mover );
        wiggle_cdr_h3->add_mover( shear_mover );
        wiggle_cdr_h3->add_mover( ccd_cycle );
        wiggle_cdr_h3->add_mover( with_flank_fold_tree );
        
        
        cdrh3_map->set_bb( flank_allow_bb_move );
        with_flank_fold_tree->apply( pose_in );
        loop_min_mover->movemap( cdrh3_map );
        loop_min_mover->apply( pose_in );
        cdrh3_map->set_bb( allow_bb_move );
        
        // rotamer trials
        select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
                             allow_repack);
        cdrh3_map->set_chi( allow_repack );
        setup_packer_task( start_pose_, tf_ );
        ( *highres_scorefxn_ )( pose_in );
        tf_->push_back( new RestrictToInterface( allow_repack ) );
        RotamerTrialsMoverOP pack_rottrial = new RotamerTrialsMover( highres_scorefxn_, tf_ );
        
        pack_rottrial->apply( pose_in );
        
        
        Real const init_temp( 2.0 );
        Real const last_temp( 0.5 );
        Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
        Real temperature = init_temp;
        
        MonteCarloOP mc;
        mc = new protocols::moves::MonteCarlo( pose_in, *highres_scorefxn_, temperature );
        mc->reset( pose_in ); // monte carlo reset
        
        bool relaxed_H3_found_ever( false );
        if( H3_filter_)
            relaxed_H3_found_ever = CDR_H3_filter( pose_in,
                                                  ab_info_->get_CDR_loop("h3")->start(),
                                                  (ab_info_->get_CDR_loop("h3")->stop()-1 - ab_info_->get_CDR_loop("h3")->start()) + 1, 
                                                  H3_filter_,
                                                  is_camelid_
                                                  );
        
        // outer cycle
        for(Size i = 1; i <= outer_cycles; i++) {
            mc->recover_low( pose_in );
            Size h3_attempts(0); // number of H3 checks after refinement
            
            // inner cycle
            for ( Size j = 1; j <= inner_cycles; j++ ) {
                temperature *= gamma;
                mc->set_temperature( temperature );
                wiggle_cdr_h3->apply( pose_in );
                cdrh3_map->set_bb( flank_allow_bb_move );
                loop_min_mover->movemap( cdrh3_map );
                loop_min_mover->apply( pose_in );
                cdrh3_map->set_bb( allow_bb_move );
                
                // rotamer trials
                select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
                                     allow_repack);
                cdrh3_map->set_chi( allow_repack );
                setup_packer_task( start_pose_, tf_ );
                ( *highres_scorefxn_ )( pose_in );
                tf_->push_back( new RestrictToInterface( allow_repack ) );
                RotamerTrialsMoverOP pack_rottrial = new RotamerTrialsMover( highres_scorefxn_, tf_ );
                pack_rottrial->apply( pose_in );
                
                bool relaxed_H3_found_current(false);
                // H3 filter check
                if(H3_filter_ && (h3_attempts <= inner_cycles)) {
                    h3_attempts++;
                    relaxed_H3_found_current = CDR_H3_filter(pose_in,
                                                             ab_info_->get_CDR_loop("h3")->start(), 
                                                             (ab_info_->get_CDR_loop("h3")->stop()-1 - ab_info_->get_CDR_loop("h3")->start()) + 1,
                                                             H3_filter_,
                                                             is_camelid_);
                    
                    if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
                        mc->boltzmann( pose_in );
                    }
                    else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
                        relaxed_H3_found_ever = true;
                        mc->reset( pose_in );
                    }
                    else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
                        --j;
                        continue;
                    }
                    else if( relaxed_H3_found_ever && relaxed_H3_found_current )
                        mc->boltzmann( pose_in );
                }
                else {
                    if( H3_filter_ ) {
                        bool relaxed_H3_found_current(false);
                        relaxed_H3_found_current = CDR_H3_filter(pose_in,
                                                                 ab_info_->get_CDR_loop("h3")->start(), 
                                                                 ( ab_info_->get_CDR_loop("h3")->stop()-1 - ab_info_->get_CDR_loop("h3")->start()) + 1,
                                                                 H3_filter_,
                                                                 is_camelid_);
                        if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
                            mc->boltzmann( pose_in );
                        }
                        else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
                            relaxed_H3_found_ever = true;
                            mc->reset( pose_in );
                        }
                        else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
                            mc->recover_low( pose_in );
                        }
                        else if( relaxed_H3_found_ever && relaxed_H3_found_current )
                            mc->boltzmann( pose_in );
                    }
                    else
                        mc->boltzmann( pose_in );
                }
                
                if ( numeric::mod(j,Size(20))==0 || j==inner_cycles ) {
                    // repack trial
                    loop_repack = new PackRotamersMover( highres_scorefxn_ );
                    setup_packer_task( start_pose_, tf_ );
                    ( *highres_scorefxn_ )( pose_in );
                    tf_->push_back( new RestrictToInterface( allow_repack ) );
                    loop_repack->task_factory( tf_ );
                    loop_repack->apply( pose_in );
                    mc->boltzmann( pose_in );
                }
            } // inner cycles
        } // outer cycles
        mc->recover_low( pose_in );
        
        // minimize
        if( !benchmark_ ) {
            cdrh3_map->set_bb( flank_allow_bb_move );
            with_flank_fold_tree->apply( pose_in );
            loop_min_mover->movemap( cdrh3_map );
            loop_min_mover->apply( pose_in );
            cdrh3_map->set_bb( allow_bb_move );
        }
        
        // Restoring pose stuff
        pose_in.fold_tree( tree_in ); // Tree
        
        TR << "H3M Finished Relaxing CDR H3 Loop" << std::endl;
        
        return;
    } //  CDRH3Modeler2::loop_fa_relax


    
    
    
    



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
 
    antibody2::Ab_Info starting_antibody;
    starting_antibody = ab_info_;
    bool closed_cutpoints( false );
 
    Size cycle( 1 );
    while( !closed_cutpoints && cycle<max_cycle_ ) {
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
 
    antibody2::Ab_Info starting_antibody;
    starting_antibody = ab_info_;
    bool closed_cutpoints( false );
 
    Size cycle (1 );
    while( !closed_cutpoints && cycle < max_cycle_ ) {
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
