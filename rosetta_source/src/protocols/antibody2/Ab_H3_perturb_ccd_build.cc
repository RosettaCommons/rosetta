// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_H3_perturb_ccd_build.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/Ab_H3_perturb_ccd_build.hh>


#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below

#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>

#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/FragmentMover.hh>


#include <protocols/antibody2/Ab_util.hh>
#include <protocols/antibody2/Ab_H3_cter_insert_mover.hh>




static numeric::random::RandomGenerator RG(21141980);




using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.Ab_H3_perturb_ccd_build");




using namespace core;
namespace protocols {
namespace antibody2 {

    
    
    
// default constructor
Ab_H3_perturb_ccd_build::Ab_H3_perturb_ccd_build() : Mover() {
    init();
}

   
    
Ab_H3_perturb_ccd_build::Ab_H3_perturb_ccd_build( bool current_loop_is_H3, bool H3_filter, bool is_camelid, Ab_Info & antibody_in ) : Mover() {
    set_default();
    init();
    
    current_loop_is_H3_ = current_loop_is_H3;
    H3_filter_ =  H3_filter;
    is_camelid_= is_camelid;
    ab_info_ = antibody_in;

}

    

    
    
void Ab_H3_perturb_ccd_build::set_default(){ 

	apply_centroid_mode_ = false;
    max_cycle_ = 20;
    c_ter_stem_ = 3;

	// size of loop above which 9mer frags are used
	cutoff_9_ = 16; // default 16

    // size of loop above which 3mer frags are used
	cutoff_3_ = 6; // default 6
    
    cen_cst_ = 10.0;
	current_loop_is_H3_ = true;
    
    H3_filter_ = true;


}
    
    
    
// default destructor
Ab_H3_perturb_ccd_build::~Ab_H3_perturb_ccd_build() {}
    
//clone
protocols::moves::MoverOP Ab_H3_perturb_ccd_build::clone() const {
    return( new Ab_H3_perturb_ccd_build() );
}
    
    

    
    
void Ab_H3_perturb_ccd_build::init( ) 
{
    if(is_camelid_) {
    }
    lowres_scorefxn_ = scoring::ScoreFunctionFactory::
    create_score_function( "cen_std", "score4L" );
	lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
	// adding constraints
	lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
    

}
    
    
void Ab_H3_perturb_ccd_build::setup_objects(){
    ab_h3_cter_insert_mover_ = new Ab_H3_cter_insert_mover(ab_info_);
    //TODO:  
    //JQX: right now, just want the code to work, whether this ab_info_ is at the right position or 
    // wehther it has been initialized or not? I don't know. Will come back to address this
}
    
    
std::string Ab_H3_perturb_ccd_build::get_name() const {
    return "Ab_H3_perturb_ccd_build";
}

    
    
    
    
void Ab_H3_perturb_ccd_build::finalize_setup( pose::Pose & pose ) {        
    read_and_store_fragments( pose );
}

    
    
    
//APPLY
void Ab_H3_perturb_ccd_build::apply( pose::Pose & pose ) {
    finalize_setup( pose );
    build_centroid_loop(pose);
}
    
    
    
    
    
    
    
    
void Ab_H3_perturb_ccd_build::build_centroid_loop( core::pose::Pose & pose ) {
        using namespace core::pose;
        using namespace core::scoring;
        using namespace protocols::moves;
        
        if( !apply_centroid_mode_ )
            return;
        
        TR <<  "H3M Modeling Centroid CDR H3 loop" << std::endl;
        
        Size frmrk_loop_end_plus_one( ab_info_.get_CDR_loop("h3")->stop() );
        Size framework_loop_size = ( frmrk_loop_end_plus_one - ab_info_.get_CDR_loop("h3")->start() ) + 1;
        Size cutpoint = ab_info_.get_CDR_loop("h3")->start() + 1;
        loops::Loop cdr_h3( ab_info_.get_CDR_loop("h3")->start(), frmrk_loop_end_plus_one,
                           cutpoint,	0, true );
        simple_one_loop_fold_tree( pose, cdr_h3 );
        
        // silly hack to make extended loops to work
        loops::LoopsOP cdr_h3_loop_list = new loops::Loops();
        cdr_h3_loop_list->add_loop( cdr_h3 );
        /* Commented out by BDW with JX's consent
         loops::loop_mover::LoopMoverOP my_loop_move = new loops::loop_mover::LoopMover( cdr_h3_loop_list );
         my_loop_move->set_extended_torsions( pose, cdr_h3 );
         my_loop_move->apply( pose );
         */
        Size unaligned_cdr_loop_begin(0), unaligned_cdr_loop_end(0);
        std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
        core::import_pose::pose_from_pdb( template_pose_, path+"hfr.pdb" );
        std::string template_name = "h3";
        antibody2::Ab_Info hfr_template( template_pose_, template_name );
        unaligned_cdr_loop_begin = hfr_template.current_start;
        unaligned_cdr_loop_end = hfr_template.current_end;
        
        pose.set_psi( ab_info_.get_CDR_loop("h3")->start() - 1,
                     template_pose_.psi( unaligned_cdr_loop_begin - 1 ) );
        pose.set_omega(ab_info_.get_CDR_loop("h3")->start() - 1,
                       template_pose_.omega( unaligned_cdr_loop_begin - 1 ) );
        
        Size modified_framework_loop_end = frmrk_loop_end_plus_one - c_ter_stem_;
        loops::Loop trimmed_cdr_h3( ab_info_.get_CDR_loop("h3")->start(),
                                   modified_framework_loop_end, cutpoint, 0, true );
        
        antibody2::Ab_Info starting_antibody;
        starting_antibody = ab_info_;
        bool closed_cutpoints( false );
        
        Size cycle ( 1 );
        while( !closed_cutpoints && cycle < max_cycle_) {
            ab_info_ = starting_antibody;
            if( framework_loop_size > 5 ){
                ab_h3_cter_insert_mover_->apply(pose);
            }
            scored_frag_close( pose, trimmed_cdr_h3 );
            if( trimmed_cdr_h3.size() > cutoff_9_  ) { // aroop_temp default cutoff_9_
                Size saved_cutoff_9 = cutoff_9_;
                cutoff_9_ = 100; // never going to reach
                scored_frag_close( pose, trimmed_cdr_h3 );
                cutoff_9_ = saved_cutoff_9; // restoring
            }
            closed_cutpoints = cutpoints_separation( pose, ab_info_ );
            ++cycle;
        } // while( ( cut_separation > 1.9 )
        
        TR <<  "H3M Finished Modeling Centroid CDR H3 loop" << std::endl;
        
        return;
        
    } // build_centroid_loop
    
    

    
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////
    /// @begin scored_frag_close
    ///
    /// @brief builds a loop from fragments file.
    ///
    /// @detailed Loop is built by a monte carlo simulation using fragments
    ///           from a fragment files. CCD moves are used to close loops
    ///           with gaps at cutpoint.H3_check is enforced if H3_filter flag
    ///           is set in command line. Loop building results in many files
    ///           containing differnt conformations of the same loop in
    ///           phi-psi-omega angles. Parallel processing is utilized.
    ///
    /// @param[in] weight_map: in this case its a centroid weight
    ///            pose_in: loop to be built on this template provided
    ///            loop_begin/loop_end: loop termini definition
    ///            frag_size: 3-mer or 9-mer
    ///            frag_offset:agreement in frag file numbering & pose numberng
    ///            cycles1: max cycles to be spent building loops
    ///            cycles2: # of fragment swaps for each loop(depends on size)
    ///            do_ccd_moves: should ccd moves be used to close gaps
    ///
    /// @global_read benchmark_
    ///
    /// @global_write
    ///
    /// @remarks
    ///
    /// @references
    ///
    /// @authors Aroop 02/04/2010
    ///
    /// @last_modified 02/04/2010
    ///////////////////////////////////////////////////////////////////////////
void Ab_H3_perturb_ccd_build::scored_frag_close (
                                           pose::Pose & pose_in,
                                           loops::Loop const trimmed_cdr_h3 ) {
        using namespace fragment;
        using namespace protocols;
        using namespace protocols::simple_moves;
        using namespace protocols::loops;
        using loop_closure::ccd::CcdMover;
        using loop_closure::ccd::CcdMoverOP;
        using loop_closure::ccd::CcdLoopClosureMover;
        using loop_closure::ccd::CcdLoopClosureMoverOP;
        
        TR <<  "H3M Fragments based centroid CDR H3 loop building" << std::endl;
        
        if( trimmed_cdr_h3.size() <= 2)
            utility_exit_with_message("Loop too small for modeling");
        
        // set cutpoint variants for correct chainbreak scoring
        if( !pose_in.residue( trimmed_cdr_h3.cut() ).is_upper_terminus() ) {
            if( !pose_in.residue( trimmed_cdr_h3.cut() ).has_variant_type(chemical::CUTPOINT_LOWER))
                core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, trimmed_cdr_h3.cut() );
            if( !pose_in.residue( trimmed_cdr_h3.cut() + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
                core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, trimmed_cdr_h3.cut() + 1 );
        }
        
        
        Size cycles1(10);
        // aroop_temp default 25 * loop size
        Size cycles2(25 * trimmed_cdr_h3.size() );
        
        // params
        Real const ccd_threshold( 0.1);
        Size h3_attempts(0);
        Real h3_fraction = 0.75; // 75% of loops are required to be H3's
        Real current_h3_prob = RG.uniform();;
        bool H3_found_ever(false);
        bool loop_found(false);
        Size total_cycles(0);
        Size frag_size(0);
        FragSetOP frags_to_use;
        {
            if( trimmed_cdr_h3.size() > cutoff_9_ ) {
                frags_to_use = cdr_h3_frags_[1]->empty_clone();
                frags_to_use = cdr_h3_frags_[1];
                frag_size = 9;
            }
            else {
                frags_to_use = cdr_h3_frags_[2]->empty_clone();
                frags_to_use = cdr_h3_frags_[2];
                frag_size = 3;
            }
        }
        
        // Storing Fold Tree
        kinematics::FoldTree old_fold_tree = pose_in.fold_tree();
        // New Fold Tree
        simple_one_loop_fold_tree( pose_in, trimmed_cdr_h3 );
        
        //setting MoveMap
        kinematics::MoveMapOP cdrh3_map;
        cdrh3_map = new kinematics::MoveMap();
        cdrh3_map->clear();
        cdrh3_map->set_chi( true );
        cdrh3_map->set_bb( false );
        for( Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop(); ii++ )
            cdrh3_map->set_bb( ii, true );
        cdrh3_map->set_jump( 1, false );
        
        // setup monte_carlo
        Real temp( 2.0);
        protocols::moves::MonteCarloOP mc, outer_mc;
        mc = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
        outer_mc = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
        Size buffer( (is_camelid_ && ab_info_.is_extended()) ? 2 : 0 );
        while( !loop_found && ( total_cycles++ < cycles1) ) {
            // insert random fragments over the whole loop
            for(Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop()
                - ( buffer + (frag_size - 1 ) ); ii++ ) {
                ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
                cfm->set_check_ss( false );
                cfm->enable_end_bias_check( false );
                cfm->define_start_window( ii );
                cfm->apply( pose_in );
            }
            if( total_cycles == 1 )
                mc->reset( pose_in );
            
            Size local_h3_attempts(0);
            for ( Size c2 = 1; c2 <= cycles2; ++c2 ) {
                // apply a random fragment
                ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
                cfm->set_check_ss( false );
                cfm->enable_end_bias_check( false );
                cfm->apply( pose_in );
                
                bool H3_found_current(false);
                if( current_loop_is_H3_ && H3_filter_ &&
                   ( local_h3_attempts++ < (50 * cycles2) ) ) {
                    H3_found_current = CDR_H3_filter(pose_in,
                                                     ab_info_.get_CDR_loop("h3")->start(),
                                                     ( ab_info_.get_CDR_loop("h3")->stop()-1 - ab_info_.get_CDR_loop("h3")->start() ) + 1,
                                                     H3_filter_,
                                                     is_camelid_);
                    if( !H3_found_ever && !H3_found_current) {
                        --c2;
                        mc->boltzmann( pose_in );
                    }
                    else if( !H3_found_ever && H3_found_current ) {
                        H3_found_ever = true;
                        mc->reset( pose_in );
                    }
                    else if( H3_found_ever && !H3_found_current ) {
                        --c2;
                        continue;
                    }
                    else if( H3_found_ever && H3_found_current )
                        mc->boltzmann( pose_in );
                }
                else
                    mc->boltzmann( pose_in );
                
                if ( (c2 > cycles2/2 && RG.uniform() * cycles2 < c2) ||
                    ( trimmed_cdr_h3.size() <= 5) ) {
                    // in 2nd half of simulation, start trying to close the loop:
                    CcdMoverOP ccd_moves = new CcdMover( trimmed_cdr_h3, cdrh3_map );
                    protocols::moves::RepeatMoverOP ccd_cycle;
                    if( trimmed_cdr_h3.size() <= 5 ) {
                        ccd_cycle = new protocols::moves::RepeatMover(ccd_moves,500*trimmed_cdr_h3.size());
                        ccd_cycle->apply( pose_in );
                    }
                    else {
                        ccd_cycle = new protocols::moves::RepeatMover(ccd_moves, 10*trimmed_cdr_h3.size());
                        ccd_cycle->apply( pose_in );
                    }
                    mc->boltzmann( pose_in );
                }
            }
            
            mc->recover_low( pose_in );
            CcdLoopClosureMoverOP ccd_closure = new CcdLoopClosureMover(
                                                                        trimmed_cdr_h3, cdrh3_map );
            ccd_closure->set_tolerance( ccd_threshold );
            ccd_closure->set_ccd_cycles( 500 );
            ccd_closure->apply( pose_in );
            
            if( total_cycles == 1 )
                outer_mc->reset( pose_in );
            
            if ( ccd_closure->forward_deviation() <= ccd_threshold &&
                ccd_closure->backward_deviation() <= ccd_threshold ) {
                // CDR-H3 filter for antibody mode
                // introduce enough diversity
                outer_mc->boltzmann( pose_in );
                if( current_loop_is_H3_ && H3_filter_ &&
                   (current_h3_prob < h3_fraction) && (h3_attempts++<50) )
                    if( !CDR_H3_filter(pose_in, 
                                       ab_info_.get_CDR_loop("h3")->start(),
                                       ( ab_info_.get_CDR_loop("h3")->stop()-1 - ab_info_.get_CDR_loop("h3")->start() ) + 1,
                                       H3_filter_,
                                       is_camelid_) 
                       )
                    {
                        continue; 
                    }
                loop_found = true;
            }
            else if( H3_filter_ )
                h3_attempts++;
        }
        outer_mc->recover_low( pose_in );
        
        // Restoring Fold Tree
        pose_in.fold_tree( old_fold_tree );
        
        TR <<  "H3M Finished Fragments based centroid CDR H3 loop building"
        << std::endl;
        
        return;
    } // scored_frag_close
    
    

    
    
    
    
void Ab_H3_perturb_ccd_build::read_and_store_fragments( core::pose::Pose & pose ) {
		using namespace chemical;
		using namespace id;
		using namespace fragment;
		using namespace core::scoring;
        

        
		// fragment initialization
		utility::vector1< FragSetOP > frag_libs;
        
		protocols::loops::read_loop_fragments( frag_libs );
        
		Size frag_size = (ab_info_.get_CDR_loop("h3")->stop()-ab_info_.get_CDR_loop("h3")->start()) + 3;
		Size cutpoint =  ab_info_.get_CDR_loop("h3")->start() + int( frag_size / 2 );
		setup_simple_fold_tree(  ab_info_.get_CDR_loop("h3")->start() - 1, cutpoint,
                               ab_info_.get_CDR_loop("h3")->stop() + 1,
                               pose.total_residue(),
                               pose );
        
		FragSetOP offset_3mer_frags;
		// a fragset of same type should be able to handle everything
		offset_3mer_frags = frag_libs[2]->empty_clone();
		FrameList loop_3mer_frames;
		Size offset = 0;
		frag_libs[2]->region_simple( 1, frag_size, loop_3mer_frames );
		for ( FrameList::const_iterator it = loop_3mer_frames.begin(),
             eit = loop_3mer_frames.end(); it!=eit; ++it ) {
			FrameOP short_frame = (*it)->clone_with_frags();
			offset++;
			short_frame->shift_to( ( ab_info_.get_CDR_loop("h3")->start() - 2 ) + offset  );
			offset_3mer_frags->add( short_frame );
		}
        
		FragSetOP offset_9mer_frags;
		// a fragset of same type should be able to handle everything
		offset_9mer_frags = frag_libs[1]->empty_clone();
		FrameList loop_9mer_frames;
		offset = 0;
		frag_libs[1]->region_simple( 1, frag_size, loop_9mer_frames );
		for ( FrameList::const_iterator it = loop_9mer_frames.begin(),
             eit = loop_9mer_frames.end(); it!=eit; ++it ) {
			FrameOP short_frame = (*it)->clone_with_frags();
			offset++;
			short_frame->shift_to( ( ab_info_.get_CDR_loop("h3")->start() - 2 ) + offset  );
			offset_9mer_frags->add( short_frame );
		}
        
		cdr_h3_frags_.push_back( offset_9mer_frags );
		cdr_h3_frags_.push_back( offset_3mer_frags );
        
		return;
} // read_and_store_fragments
    

    
    
    
    
    
    
    
    
    
    
    
    
}
}



