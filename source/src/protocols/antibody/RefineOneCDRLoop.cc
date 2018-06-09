// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/RefineOneCDRLoop.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#include <protocols/antibody/RefineOneCDRLoop.hh>
#include <protocols/antibody/RefineOneCDRLoopCentroid.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <core/pose/Pose.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/H3RefineCCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>

#include <core/chemical/VariantType.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.antibody.RefineOneCDRLoop" );

using namespace core;
namespace protocols {
namespace antibody {

// default constructor
RefineOneCDRLoop::RefineOneCDRLoop( ) : Mover() {}

RefineOneCDRLoop::RefineOneCDRLoop( AntibodyInfoOP antibody_info ) : Mover() {
	user_defined_ = false;
	preset_foldtree_ = false;
	ab_info_ = antibody_info;
	init();
}

RefineOneCDRLoop::RefineOneCDRLoop( loops::Loop loop, std::string refine_mode, scoring::ScoreFunctionCOP scorefxn ) : Mover() {
	// warning if you are directly passing a loop object, I ASSUME YOUR FOLD TREE IS SET PROPERLY
	user_defined_ = true;
	one_cdr_loop_ = loop;
	preset_foldtree_ = true;
	refine_mode_  = refine_mode;
	scorefxn_ = scorefxn->clone();
	init();
}

RefineOneCDRLoop::RefineOneCDRLoop( AntibodyInfoOP antibody_info, std::string refine_mode ) : Mover() {
	user_defined_ = false;
	preset_foldtree_ = false;
	ab_info_    = antibody_info;
	refine_mode_  = refine_mode;
	init();
}


RefineOneCDRLoop::RefineOneCDRLoop(AntibodyInfoOP antibody_info,
	std::string refine_mode,
	scoring::ScoreFunctionCOP scorefxn
) : Mover()  {
	user_defined_ = true;
	preset_foldtree_ = false;
	cdr_loop_name_ = h3;
	ab_info_    = antibody_info;
	one_cdr_loop_ = ab_info_->get_CDR_loop(cdr_loop_name_);
	refine_mode_  = refine_mode;
	scorefxn_ = scorefxn->clone();

	init();
}


RefineOneCDRLoop::RefineOneCDRLoop(AntibodyInfoOP antibody_info,
	CDRNameEnum const & cdr_loop_name,
	std::string refine_mode,
	scoring::ScoreFunctionCOP scorefxn
) : Mover()  {
	user_defined_ = true;
	preset_foldtree_ = false;
	ab_info_    = antibody_info;
	cdr_loop_name_ = cdr_loop_name;
	one_cdr_loop_ = ab_info_->get_CDR_loop(cdr_loop_name_);
	refine_mode_  = refine_mode;
	scorefxn_ = scorefxn->clone();

	init();
}


void RefineOneCDRLoop::init( ) {
	set_default();
}


void RefineOneCDRLoop::set_default() {
	flank_size_     = 2;
	flank_relax_    = true;
	high_cst_       = 100.0;
	H3_filter_      = false;
	num_filter_tries_ = 20;

	if ( !user_defined_ ) {
		refine_mode_    = "legacy_refine_ccd";
		cdr_loop_name_  = h3;
		one_cdr_loop_ = ab_info_->get_CDR_loop(cdr_loop_name_);
		scorefxn_ = scoring::get_score_function();
		scorefxn_->set_weight( scoring::chainbreak, 1.0 );
		scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
		scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );
	}
}


// default destructor
RefineOneCDRLoop::~RefineOneCDRLoop() = default;

//clone
protocols::moves::MoverOP RefineOneCDRLoop::clone() const {
	return( protocols::moves::MoverOP( new RefineOneCDRLoop() ) );
}


std::string RefineOneCDRLoop::get_name() const {
	return "RefineOneCDRLoop";
}

void RefineOneCDRLoop::set_score_function(core::scoring::ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn->clone();
}


void RefineOneCDRLoop::pass_start_pose(core::pose::Pose & start_pose) {
	start_pose_ = start_pose;
}


void RefineOneCDRLoop::apply(core::pose::Pose &pose) {

	//JQX: make sure the pose and the scoring function is on the same resolution
	if ( pose.is_fullatom() != scorefxn_->has_nonzero_weight( core::scoring::fa_rep ) ) {
		utility_exit_with_message("the resultions of the 'pose' and the 'scoring function' don't match!");
	}

	// add cutpoint variants (necessary for chainbreak score)
	loops::remove_cutpoint_variants( pose, true );
	loops::add_cutpoint_variants( pose );

	if ( refine_mode_ == "legacy_centroid_refine_ccd" ) {
		if ( !preset_foldtree_ ) {
			RefineOneCDRLoopCentroidOP legacy_centroid_refine_ccd( new RefineOneCDRLoopCentroid( ab_info_, cdr_loop_name_, scorefxn_ ) );
			legacy_centroid_refine_ccd -> apply(pose);
		} else {
			utility_exit_with_message("In RefineOneCDRLoop, \"legacy_centroid_refine_ccd\" option is in compatible with \"preset_foldtree_\"!");
		}
	} else if ( refine_mode_ == "legacy_refine_ccd" ) {  //the legacy centroid_refine_ccd method
		if ( !preset_foldtree_ ) {
			H3RefineCCDOP legacy_refine_ccd( new H3RefineCCD(ab_info_, cdr_loop_name_, scorefxn_) );
			legacy_refine_ccd -> pass_start_pose(start_pose_);
			if ( !flank_relax_ ) {
				legacy_refine_ccd->turn_off_flank_relax();
			}
			legacy_refine_ccd -> apply(pose);
		} else {
			utility_exit_with_message("In RefineOneCDRLoop, \"legacy_refine_ccd\" option is in compatible with \"preset_foldtree_\"!");
		}
	} else {  // options do not entail legacy methods ...
		if ( !preset_foldtree_ ) {
			// do some work on the fold tree
			// note this is different-ish from the loop we passed
			loops::Loop one_cdr_loop = ab_info_->get_CDR_loop(cdr_loop_name_);

			if ( flank_relax_ ) {
				// JQX: The idea is to minimize the flanking residues backbones on each stems,
				//      but only minimization, no perturbation (small, shear, etc ..) to them.
				//      The fold_tree should be set correctly for the MinMover (or AtomTreeMinimizer).
				//      Will this special fold_tree affect perturbation? I don't think so. It doesn't really matter
				//      2 or 3 residues before loop_begin or after loop_end, as long as the cut point is in
				//      the middle and the loop is within two jumps points.
				simple_fold_tree( pose,  one_cdr_loop.start()- 1 - flank_size_, one_cdr_loop.cut(), one_cdr_loop.stop() + 1 + flank_size_ );
			} else {
				simple_fold_tree( pose,  one_cdr_loop.start()- 1, one_cdr_loop.cut(), one_cdr_loop.stop() + 1 );
			}
			//JQX: be careful, no fold_tree operations are conducted in the two movers below.
		}

		// LoopMovers take LoopsOp so wrap one CDR loop
		loops::LoopsOP pass_loops( new loops::Loops() );
		pass_loops->add_loop(  one_cdr_loop_  );


		core::Size itry=1;
		core::Real best_score=0.0;
		core::pose::Pose best_pose;
		std::stringstream Num;
		std::string str_num;
		core::pose::Pose pose_before_refine = pose; //JQX: save the current pose before doing any refinement

		if ( refine_mode_ == "refine_ccd" ) {

			loops::loop_mover::refine::LoopMover_Refine_CCD refine_ccd( pass_loops, scorefxn_ );

			if ( get_native_pose() ) {
				refine_ccd.set_native_pose( get_native_pose() );
			}

			if ( flank_relax_ ) {
				refine_ccd.set_flank_residue_min(flank_relax_);
			}

			while ( itry<=num_filter_tries_ ) {
				TR<<"   Trying Refinement  ................. "<<itry<<std::endl;
				pose = pose_before_refine;
				refine_ccd.apply( pose );

				Num<<itry;
				Num>>str_num;
				Num.str("");
				Num.clear();

				if ( H3_filter_ ) { // this is deprecated and should not be true
					if ( CDR_H3_cter_filter(pose, ab_info_) ) {
						break;
					} else {
						scorefxn_->score(pose); //Segfault protection
						if (  pose.energies().total_energy() <= best_score ) {
							best_score = pose.energies().total_energy();
							best_pose = pose;
						}
					}
				} else { // H3_filter_ == false
					break;
				}
				if ( itry==num_filter_tries_ ) {
					pose=best_pose;
				}
				itry++;
			}

		} else if ( refine_mode_ == "refine_kic" ) {

			loops::loop_mover::refine::LoopMover_Refine_KIC refine_kic( pass_loops, scorefxn_ );

			if ( get_native_pose() ) {
				refine_kic.set_native_pose( get_native_pose() );
			}

			if ( flank_relax_ )      {
				refine_kic.set_flank_residue_min(flank_relax_);
			}

			while ( itry<=num_filter_tries_ ) {
				TR<<"   Trying Refinement  ................. "<<itry<<std::endl;
				pose = pose_before_refine;
				refine_kic.apply( pose );

				Num<<itry;
				Num>>str_num;
				Num.str("");
				Num.clear();

				if ( H3_filter_ ) { // this is deprecated and should not be true
					if ( CDR_H3_cter_filter(pose, ab_info_) ) {
						break;
					} else {
						scorefxn_->score(pose); //Segfault protection
						if (  pose.energies().total_energy() <= best_score ) {
							best_score = pose.energies().total_energy();
							best_pose = pose;
						}
					}
				} else { // H3_filter_ == false
					break;
				}
				if ( itry==num_filter_tries_ ) {
					pose=best_pose;
				}
				itry++;
			}

		} else {
			utility_exit_with_message("the refinement method is not available!");
		}

	}// refine_CCD or refine_KIC


}//apply


} // namespace antibody
} // namespace protocols


