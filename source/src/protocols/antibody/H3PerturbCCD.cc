// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody/H3PerturbCCD.cc
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody/H3PerturbCCD.hh>

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

#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>

#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>



static numeric::random::RandomGenerator RG(21141980);




using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody.H3PerturbCCD");




using namespace core;
namespace protocols {
namespace antibody {




// default constructor
H3PerturbCCD::H3PerturbCCD() : Mover() {
	user_defined_ = false;
}

// default destructor
H3PerturbCCD::~H3PerturbCCD() {}


H3PerturbCCD::H3PerturbCCD( AntibodyInfoOP antibody_in ) : Mover() {
	user_defined_ = false;
	ab_info_ = antibody_in;

	init();
}


H3PerturbCCD::H3PerturbCCD( AntibodyInfoOP antibody_in,
                            core::scoring::ScoreFunctionCOP lowres_scorefxn ) : Mover() {
	user_defined_ = true;
	ab_info_=antibody_in;
	lowres_scorefxn_ = lowres_scorefxn->clone();

	init();
}



void H3PerturbCCD::init( ) {
	set_default();

}




void H3PerturbCCD::set_default() {
	cutoff_9_       = 16; // size of loop above which 9mer frags are used
	cutoff_3_       = 6;  // size of loop above which 3mer frags are used
	cen_cst_        = 10.0;
	num_cycles1_    = 10;  //max cycles to be spent building loops
	h3_fraction_    = 0.75; // 75% of loops are required to be H3's
	max_ccd_cycles_ = 500;
	ccd_threshold_  = 0.1;
	Temperature_    = 2.0;

	is_camelid_         = false;
	current_loop_is_H3_ = true;
	H3_filter_          = true;



	if(!user_defined_) {
		lowres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
		lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
		lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
	}
}





//clone
protocols::moves::MoverOP H3PerturbCCD::clone() const {
	return( new H3PerturbCCD() );
}







void H3PerturbCCD::finalize_setup( pose::Pose & pose_in ) {

	read_and_store_fragments(  );

	mc_ = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, Temperature_ );
	outer_mc_ = new protocols::moves::MonteCarlo( pose_in, *lowres_scorefxn_, Temperature_ );

}




void H3PerturbCCD::apply( pose::Pose & pose_in ) {


	finalize_setup( pose_in );


	loops::Loop trimmed_cdr_h3 = input_loop_;

	using namespace fragment;
	using namespace protocols;
	using namespace protocols::simple_moves;
	using namespace protocols::loops;
	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	TR <<  "Fragments based centroid CDR H3 loop building" << std::endl;

	if( trimmed_cdr_h3.size() <= 2) {
		utility_exit_with_message("Loop too small for modeling");
	}

	// params
	Size h3_attempts(0);
	Real current_h3_prob = RG.uniform();
	TR<<"current_h3_prob="<<current_h3_prob<<std::endl;

	Size frag_size(0);
	FragSetOP frags_to_use;
	{
		if( trimmed_cdr_h3.size() > cutoff_9_ ) {
			frags_to_use = cdr_h3_frags_[1]->empty_clone();
			frags_to_use = cdr_h3_frags_[1];
			frag_size = 9;
		} else {
			frags_to_use = cdr_h3_frags_[2]->empty_clone();
			frags_to_use = cdr_h3_frags_[2];
			frag_size = 3;
		}
	}
	TR<<"frag_size="<<frag_size<<std::endl;

	// Storing Fold Tree
	kinematics::FoldTree old_fold_tree = pose_in.fold_tree();

	// New Fold Tree
	simple_one_loop_fold_tree( pose_in, trimmed_cdr_h3 );// is the cutpoint important or not???
	TR<<trimmed_cdr_h3<<std::endl;
	TR<<pose_in<<std::endl;


	// set cutpoint variants for correct chainbreak scoring
	loops::remove_cutpoint_variants( pose_in, true );
	if( !pose_in.residue( trimmed_cdr_h3.cut() ).is_upper_terminus() ) {
		if( !pose_in.residue( trimmed_cdr_h3.cut() ).has_variant_type(chemical::CUTPOINT_LOWER))
			core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_LOWER, trimmed_cdr_h3.cut() );
		if( !pose_in.residue( trimmed_cdr_h3.cut() + 1 ).has_variant_type(chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose_in, chemical::CUTPOINT_UPPER, trimmed_cdr_h3.cut() + 1 );
	}


	//setting MoveMap
	//JQX: all the chi angles of all the side chains are flexible
	//     only the backbone of the trimmed_cdr_h3 is flexible
	kinematics::MoveMapOP cdrh3_map;
	cdrh3_map = new kinematics::MoveMap();
	cdrh3_map->clear();
	cdrh3_map->set_chi(true );
	cdrh3_map->set_bb (false);
	for( Size ii=trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop(); ii++ ) {
		cdrh3_map->set_bb( ii, true );
	}
	cdrh3_map->set_jump( 1, false );


	// aroop_temp default 25 * loop size
	Size num_cycles2(25 * trimmed_cdr_h3.size() );
	bool H3_found_ever(false);
	Size total_cycles(0);
	Size buffer(   (is_camelid_ && (ab_info_->get_H3_kink_type()==Extended)   ) ? 2 : 0 );
	bool loop_found(false);

	while( !loop_found && ( total_cycles++ < num_cycles1_) ) {
		// JQX: insert random fragments over the whole loop
		//      fragment has a size frag_size
		//TR<<"trimmed_cdr_h3.start() = "<<trimmed_cdr_h3.start()<<std::endl;
		//TR<<"trimmed_cdr_h3.stop() - ( buffer + (frag_size - 1 ) ) = "<<trimmed_cdr_h3.stop() - ( buffer + (frag_size - 1 ) )<<std::endl;
		for(Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop() - ( buffer + (frag_size - 1 ) ); ii++ ) {
			ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			cfm->define_start_window( ii );
			cfm->apply( pose_in );
		}


		if( total_cycles == 1 ) {
			mc_->reset( pose_in );
		}
		Size local_h3_attempts(0);

		for ( Size c2 = 1; c2 <= num_cycles2; ++c2 ) {
			TR<<"c1="<<total_cycles<<"    "<<"c2="<<c2<<std::endl;
			// apply a random fragment
			ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use, cdrh3_map);
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			cfm->apply( pose_in );

			bool H3_found_current(false);
			if( current_loop_is_H3_ && H3_filter_ && ( local_h3_attempts++ < (50 * num_cycles2) ) ) {
				H3_found_current = CDR_H3_cter_filter(pose_in, ab_info_);
				if( !H3_found_ever && !H3_found_current) {
					--c2;
					mc_->boltzmann( pose_in );
					// JQX: 1. this try failed, and never succeeded before, --c2
					//      2. accept or reject this pose, then retry this cycle
				} else if( !H3_found_ever && H3_found_current ) {
					H3_found_ever = true;
					mc_->reset( pose_in );
					// JQX: 1. this try succeeded, and it is the 1st time success
					//      2. reset the pose
				} else if( H3_found_ever && !H3_found_current ) {
					--c2;
					continue;
					// JQX: this try failed, but you get success before, --c2, retry this cycle
				} else if( H3_found_ever && H3_found_current ) {
					mc_->boltzmann( pose_in );
					// JQX: this try succeeded, and you also succeeded before
				}
			} else {
				mc_->boltzmann( pose_in );
			}

			// TODO:
			// JQX: this "RG.uniform() * num_cycles2 < c2" is so weird, not sure what Aroop really wants to do
			if ( (c2 > num_cycles2/2 && RG.uniform() * num_cycles2 < c2) || ( trimmed_cdr_h3.size() <= 5) ) {
				// in 2nd half of simulation, start trying to close the loop:
				CCDLoopClosureMoverOP ccd_moves = new CCDLoopClosureMover( trimmed_cdr_h3, cdrh3_map );
				protocols::moves::RepeatMoverOP ccd_cycle;
				if( trimmed_cdr_h3.size() <= 5 ) {
					ccd_cycle = new protocols::moves::RepeatMover(ccd_moves,500*trimmed_cdr_h3.size());
					ccd_cycle->apply( pose_in );
				} else {
					ccd_cycle = new protocols::moves::RepeatMover(ccd_moves, 10*trimmed_cdr_h3.size());
					ccd_cycle->apply( pose_in );
				}
				mc_->boltzmann( pose_in );
			}
		}// finish cycles2


		mc_->recover_low( pose_in );

		CCDLoopClosureMoverOP ccd_closure = new CCDLoopClosureMover(trimmed_cdr_h3, cdrh3_map );
		ccd_closure->tolerance( ccd_threshold_ );
		ccd_closure->max_cycles( max_ccd_cycles_ );
		ccd_closure->apply( pose_in );

		if( total_cycles == 1 ) {
			outer_mc_->reset( pose_in );
		}

		if ( ccd_closure->deviation() <= ccd_threshold_ ) {
			// CDR-H3 filter for antibody mode
			// introduce enough diversity
			outer_mc_->boltzmann( pose_in );
			if( current_loop_is_H3_ && H3_filter_ && (current_h3_prob < h3_fraction_) && (h3_attempts++<50) ) {
				if(   !CDR_H3_cter_filter(pose_in, ab_info_)    ) {
					continue;
				}
			}
			loop_found = true;
		} else if( H3_filter_ ) {
			h3_attempts++;
		}

	}// finish cycles1

	outer_mc_->recover_low( pose_in );

	// Restoring Fold Tree
	pose_in.fold_tree( old_fold_tree );

	TR <<  "Finished Fragments based centroid CDR H3 loop building" << std::endl;

	return;
}






void H3PerturbCCD::read_and_store_fragments( ) {
	using namespace chemical;
	using namespace id;
	using namespace fragment;
	using namespace core::scoring;



	// fragment initialization
	utility::vector1< FragSetOP > frag_libs;

	protocols::loops::read_loop_fragments( frag_libs );

	Size frag_size = (ab_info_->get_CDR_loop(h3).stop() - ab_info_->get_CDR_loop(h3).start()) + 3; //JQX: why +3??
	TR<<frag_size<<std::endl;



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
		short_frame->shift_to( ( ab_info_->get_CDR_loop(h3).start() - 2 ) + offset  );
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
		short_frame->shift_to( ( ab_info_->get_CDR_loop(h3).start() - 2 ) + offset  );
		offset_9mer_frags->add( short_frame );
	}



	cdr_h3_frags_.push_back( offset_9mer_frags );
	cdr_h3_frags_.push_back( offset_3mer_frags );

	TR<<"Finished reading fragments files!!!"<<std::endl;

	return;

}










}// namespace antibody
}// namespace protocols



