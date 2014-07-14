// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/IndependentLoopMover.cc
/// @brief  loop mover base class
/// @author Mike Tyka
/// @author James Thompson

/// Unit headers
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/loops_main.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/checkpoint/CheckPointer.hh>

// Basic headers
#include <basic/Tracer.hh> // tracer output
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random_permutation.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


namespace protocols {
namespace loops {
namespace loop_mover {

static numeric::random::RandomGenerator RG(47537); // <- Magic number, do not change it!!!

///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

IndependentLoopMover::IndependentLoopMover() :
	LoopMover()
{
	Mover::type("IndependentLoopMover");
	set_defaults();
}

IndependentLoopMover::IndependentLoopMover(
	protocols::loops::LoopsOP loops_in
) : LoopMover( loops_in )
{
	Mover::type("IndependentLoopMover");
	set_defaults();
}

IndependentLoopMover::IndependentLoopMover(
  utility::vector1< bool > const& selection
) : LoopMover( new protocols::loops::Loops( selection ) )
{
  Mover::type("IndependentLoopMover");
  set_defaults();
}

IndependentLoopMover::IndependentLoopMover(
	protocols::loops::LoopsFileData const & lfd
) : LoopMover( lfd )
{
	Mover::type("IndependentLoopMover");
	set_defaults();
}

IndependentLoopMover::IndependentLoopMover(
	GuardedLoopsFromFileOP guarded_loops
) : LoopMover( guarded_loops )
{
	Mover::type("IndependentLoopMover");
	set_defaults();
}

// destructor
IndependentLoopMover::~IndependentLoopMover(){}

void IndependentLoopMover::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	build_attempts_        = option[ OptionKeys::loops::build_attempts ]();  //     3
	grow_attempts_         = option[ OptionKeys::loops::grow_attempts ]();  // 7
	accept_aborted_loops_  = option[ OptionKeys::loops::accept_aborted_loops]();  //false
	strict_loops_          = option[ OptionKeys::loops::strict_loops ]();  //false
	random_order_          = option[ OptionKeys::loops::random_order ]();     //false
	build_all_loops_       = option[ OptionKeys::loops::build_all_loops ]();  //false
	loop_combine_rate_     = option[ OptionKeys::loops::combine_rate ]();     // 0.0
}

/// @brief Apply the loop-build protocol to the input pose
void IndependentLoopMover::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	resolve_loop_indices( pose );

 	// Select Loops to be built
	all_loops_closed_ = true;
	tr().Info << "ALL_LOOPS:" << *loops() << std::endl;

	Loops selected_loops;
	select_loops( selected_loops );
	tr().Info << "SELECTEDLOOPS:" << selected_loops << std::endl;


	kinematics::FoldTree f_orig=pose.fold_tree();
	Size lcount=0;

	core::pose::Pose pose_initial = pose;


	int select_best_loop_from = option[ OptionKeys::loops::select_best_loop_from ]();

	LoopResult result_of_loopModel;

	for ( Loops::iterator it=selected_loops.v_begin(), it_end=selected_loops.v_end();
		 it != it_end; ++it ) {
		lcount++;
		// Make loal copy of loop to be build
		Loop buildloop( *it );

		// either extend or at least idealize the loop (just in case).
		if ( buildloop.is_extended() ){
			// store starting fold tree and cut pose_initial
			set_single_loop_fold_tree( pose_initial, buildloop );
			tr().Info << "Setting extended torsions: " << buildloop << std::endl;
			if (  option[ OptionKeys::loops::debug ]() ) pose_initial.dump_pdb("just_before_set_extended_torsions.pdb");
			set_extended_torsions( pose_initial, buildloop );
			if (  option[ OptionKeys::loops::debug ]() ) pose_initial.dump_pdb("just_after_set_extended_torsions.pdb");
			pose_initial.fold_tree( f_orig );
		}

		// statistics:
		int  time_start = time(NULL);
		Size nfailure = 0;
		//Size nrmsfail = 0;

		tr().Info << "Building Loop: " << buildloop << std::endl;
		result_of_loopModel = Failure;

		pose::Pose best_pose = pose_initial;
		Real       best_score = 10000000.0;
		Size       best_count = 0;
		// code below goes for grow_attempts + 1 ...
		for ( int extension_attempt = 0; extension_attempt <= grow_attempts_; extension_attempt ++ ){

			if ( !strict_loops_ && extension_attempt > 0 ){
                loops()->grow_loop_away_from_sheets( pose, buildloop, 1.0 );
			}
			for ( int build_attempt = 0; build_attempt < build_attempts_; build_attempt ++ ){
				tr().Info << "Building Loop attempt: " << build_attempt << std::endl;
				pose = pose_initial;

				if (  option[ OptionKeys::loops::debug ]() ) pose.dump_pdb("just_before_rebuild.pdb");

				std::string checkname = "loop_" + string_of( lcount ) + "_" + string_of( extension_attempt ) + "_" + string_of( build_attempt );

				std::string curr_job_tag = get_current_tag();
				bool checkpoint_recovery = false;

				if ( get_checkpoints()->recover_checkpoint( pose, curr_job_tag, checkname + "_S", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result_of_loopModel = Success;
				} else if ( get_checkpoints()->recover_checkpoint( pose, curr_job_tag, checkname + "_C", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result_of_loopModel = CriticalFailure;
				} else if ( get_checkpoints()->recover_checkpoint( pose, curr_job_tag, checkname + "_F", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result_of_loopModel = Failure;
				} else {
					// this should have been called before here, but there are
					// some cases where loop-building is attempted with a cut
					// of zero.
					buildloop.auto_choose_cutpoint( pose );
					/// Main loop modeling call here.  Derived classes override the "model_loop" function"
					result_of_loopModel = model_loop( pose, buildloop );
				}

				// If the loop bas built and closed ok.
				if ( result_of_loopModel == Success || accept_aborted_loops_ ){
					if ( ! checkpoint_recovery ){
						get_checkpoints()->checkpoint( pose, curr_job_tag, checkname + "_S", true );
					}

					// briefly score the pose with whatever scorefunction is being used.
					core::Real pose_score;
					if ( scorefxn() ) pose_score = (*scorefxn())(pose);
					else              pose_score = 0;

					get_checkpoints()->debug(  curr_job_tag, checkname, pose_score);

					// compare to previous score - if better "accept"
					if ( pose_score < best_score || best_count == 0 ){
						best_pose = pose;
						best_score = pose_score;
						best_count ++;
						tr().Debug << "Adding a " << best_score << std::endl;
					}
					if ( best_count >= (Size)select_best_loop_from ) break;
					continue;
				}
				nfailure++;

				//fpd if we have strict loops on, keep trying rather than immediately failing
				if ( result_of_loopModel == ExtendFailure && !strict_loops_ ){ // means extend loop immediately!
					if ( ! checkpoint_recovery ){
						get_checkpoints()->checkpoint( pose, curr_job_tag, checkname + "_F", true );
					}
					get_checkpoints()->debug(  curr_job_tag, checkname, -1);
					break;
				}

				if ( result_of_loopModel == CriticalFailure ){
					if ( ! checkpoint_recovery ){
						get_checkpoints()->checkpoint( pose, curr_job_tag, checkname + "_C", true );
					}
					get_checkpoints()->debug(  curr_job_tag, checkname, -2);
					tr().Error << "Unable to build this loop - a critical error occured. Moving on .. " << std::endl;
					break;
				}
			} // for build_attempts

			// If we have a sufficient number of loop successses, break (by
			// default this means 1 closed loop)
			if ( best_count >= (Size)select_best_loop_from ) break;
			// If we can't build this loop due to some major fundamental failure
			if ( result_of_loopModel == CriticalFailure ) break;
			// or if we're still unsuccessful (i.e. best_count == 0) and growing the loop isnt an option (because strict_loops is set) then also give up.
			if ( strict_loops_ ) break;
		}

		tr().Info << "result of loop closure:0 success, 3 failure " << result_of_loopModel << std::endl;
		if(result_of_loopModel != Success){
			all_loops_closed_ = false;
			break; //no need to check the rest of the loops if one can't be closed the result will be an open structure.
		}

		if ( (best_count > 0) || accept_aborted_loops_ ){
			pose_initial = best_pose;
			pose = best_pose;
		}


		// Print statistics:
		int time_end = time(NULL);
		float time_per_build = float(time_end - time_start) / float(nfailure+1);

		tr().Info   << "Loopstat: "
			<< "  " << I(3,it->start())
			<< "  " << I(3,it->stop())
			<< "  " << I(3,buildloop.start() )
			<< "  " << I(3,buildloop.stop() )
			<< "  " << I(3,buildloop.size())
			<< "  time " << F(5,1,time_per_build )
			<< "  " << I(3,nfailure)
			<< "  time " << F(5,1,time_per_build * nfailure )
			<< "  " << best_score
			<< "  " << ( it->is_extended() ? std::string(" ext ") : std::string(" noext "))  << std::endl;

	}

	loops::remove_cutpoint_variants( pose );
	pose.fold_tree( f_orig );
}

void IndependentLoopMover::select_loops( Loops & selected_loops ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Loops temp_loops;
	Loops comb_loops;

	if ( option[ OptionKeys::loops::build_specific_loops ].user() ) {
		// Choose loops by user
		utility::vector1<int> loop_numbers(
			option[ OptionKeys::loops::build_specific_loops ]
		);

		for ( Size i = 1; i <= loop_numbers.size(); ++i ) {
				if ( loop_numbers[i]  <= 0 ){
					utility_exit_with_message( "Specified loop numbers is 0 or below.");
				}
				if ( (Size)loop_numbers[i]  > loops()->size()  ){
					utility_exit_with_message( "Specified loop number is greater than the number of loops in the loop file itself.");
				}
				temp_loops.add_loop( ( *loops() )[ loop_numbers[i] ] );
		}

	} else {
		// Choose loops by skiprate
		for ( Loops::const_iterator it=loops()->begin(), it_end=loops()->end(); it != it_end; ++it ) {
			if ( build_all_loops_ ||
				( numeric::random::uniform() >= it->skip_rate() )
			) {
				temp_loops.add_loop( *it );
			}
		}
	}

  // combine loops sing combine_rate
	std::sort( temp_loops.v_begin(), temp_loops.v_end(), Loop_lt() );  // necessary

  for ( Size l=1; l <= temp_loops.size(); l++  ) {
		if ( ( l == temp_loops.size() ) ||
			( numeric::random::uniform() >= loop_combine_rate_ )
		) {
			comb_loops.add_loop( temp_loops[ l ] );
		} else {
			Loop combined_loop(
				temp_loops[ l ].start(),
				temp_loops[ l+1 ].stop(),
				temp_loops[ l ].cut(),
				temp_loops[ l ].skip_rate() * temp_loops[ l+1 ].skip_rate(),
				temp_loops[ l ].is_extended() || temp_loops[ l+1 ].is_extended()
			);

			comb_loops.add_loop( combined_loop );
			l++; // skip next loop;
		}
	}

	// randomize order if required
	if ( random_order_ ) {
		//std::random__shuffle( comb_loops.v_begin(), comb_loops.v_end() );
		numeric::random::random_permutation( comb_loops.v_begin(), comb_loops.v_end(), RG);
	}

	selected_loops = comb_loops;
} // IndependentLoopMover::select_loops()

std::string
IndependentLoopMover::get_name() const {
	return "IndependentLoopMover";
}

} // namepsace loop_mover
} // namespace loops
} // namespace protocols
