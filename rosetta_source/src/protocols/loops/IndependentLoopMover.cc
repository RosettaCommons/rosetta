// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/IndependentLoopMover.cc
/// @brief  loop mover base class
/// @author Mike Tyka
/// @author James Thompson

#include <protocols/loops/IndependentLoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh> // tracer output

#include <core/fragment/FragSet.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/checkpoint/CheckPointer.hh>

//Utility Headers
// AUTO-REMOVED #include <numeric/random/random.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif
// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/types.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto Headers





namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

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
	basic::Tracer tr("protocol.loops.LoopMover");

 	// Select Loops to be built

	tr.Info << "ALL_LOOPS:" << loops_ << std::endl;

	Loops selected_loops;
	select_loops( selected_loops );
	tr.Info << "SELECTEDLOOPS:" << selected_loops << std::endl;


	kinematics::FoldTree f_orig=pose.fold_tree();
	Size lcount=0;

	core::pose::Pose pose_initial = pose;


	int select_best_loop_from = option[ OptionKeys::loops::select_best_loop_from ]();

	for ( Loops::iterator it=selected_loops.v_begin(), it_end=selected_loops.v_end();
		 it != it_end; ++it ) {
		lcount++;
		// Make loal copy of loop to be build
		Loop buildloop( *it );

		// either extend or at least idealize the loop (just in case).
		if ( buildloop.is_extended() ){
			// store starting fold tree and cut pose_initial
			set_single_loop_fold_tree( pose_initial, buildloop );
			tr.Info << "Setting extended torsions: " << buildloop << std::endl;
			if (  option[ OptionKeys::loops::debug ]() ) pose_initial.dump_pdb("just_before_set_extended_torsions.pdb");
			set_extended_torsions( pose_initial, buildloop );
			if (  option[ OptionKeys::loops::debug ]() ) pose_initial.dump_pdb("just_after_set_extended_torsions.pdb");
			pose_initial.fold_tree( f_orig );
		}

		// statistics:
		int  time_start = time(NULL);
		Size nfailure = 0;
		//Size nrmsfail = 0;

		tr.Info << "Building Loop: " << buildloop << std::endl;
		LoopResult result = Failure;

		pose::Pose best_pose = pose_initial;
		Real       best_score = 10000000.0;
		Size       best_count = 0;
		// code below goes for grow_attempts + 1 ...
		for ( int extension_attempt = 0; extension_attempt <= grow_attempts_; extension_attempt ++ ){

			if ( !strict_loops_ && extension_attempt > 0 ){
				loops_.grow_loop( pose, buildloop, 1.0 );
			}
			for ( int build_attempt = 0; build_attempt < build_attempts_; build_attempt ++ ){
				tr.Info << "Building Loop attempt: " << build_attempt << std::endl;
				pose = pose_initial;

				if (  option[ OptionKeys::loops::debug ]() ) pose.dump_pdb("just_before_rebuild.pdb");

				std::string checkname = "loop_" + string_of( lcount ) + "_" + string_of( extension_attempt ) + "_" + string_of( build_attempt );

				std::string curr_job_tag = get_current_tag();
				bool checkpoint_recovery = false;

				if ( checkpoints_.recover_checkpoint( pose, curr_job_tag, checkname + "_S", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result = Success;
				}
				else
				if ( checkpoints_.recover_checkpoint( pose, curr_job_tag, checkname + "_C", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result = CriticalFailure;
				}else
				if ( checkpoints_.recover_checkpoint( pose, curr_job_tag, checkname + "_F", pose.is_fullatom(), true) ) {
					checkpoint_recovery = true;
					result = Failure;
				} else {
					// this should have been called before here, but there are
					// some cases where loop-building is attempted with a cut
					// of zero.
					buildloop.auto_choose_cutpoint( pose );
					result = model_loop( pose, buildloop );
				}

				// If the loop bas built and closed ok.
				if ( result == Success || accept_aborted_loops_ ){
					if ( ! checkpoint_recovery ){
						checkpoints_.checkpoint( pose, curr_job_tag, checkname + "_S", true );
					}
					// briefly score the pose with whatever scorefunction is being used.
					core::Real pose_score;
					if ( scorefxn_ ) pose_score = (*scorefxn_)(pose);
					else             pose_score = 0;

					checkpoints_.debug(  curr_job_tag, checkname, pose_score);

					// compare to previous score - if better "accept"
					if ( pose_score < best_score || best_count == 0 ){
						best_pose = pose;
						best_score = pose_score;
						best_count ++;
						tr.Debug << "Adding a " << best_score << std::endl;
					}
					if ( best_count >= (Size)select_best_loop_from ) break;
					continue;
				}
				nfailure++;

				//fpd if we have strict loops on, keep trying rather than immediately failing
				if ( result == ExtendFailure && !strict_loops_ ){ // means extend loop immediately!
					if ( ! checkpoint_recovery ){
						checkpoints_.checkpoint( pose, curr_job_tag, checkname + "_F", true );
					}
					checkpoints_.debug(  curr_job_tag, checkname, -1);
					break;
				}

				if ( result == CriticalFailure ){
					if ( ! checkpoint_recovery ){
						checkpoints_.checkpoint( pose, curr_job_tag, checkname + "_C", true );
					}
					checkpoints_.debug(  curr_job_tag, checkname, -2);
					tr.Error << "Unable to build this loop - a critical error occured. Moving on .. " << std::endl;
					break;
				}
			} // for build_attempts

			// If we have a sufficient number of loop successses, break (by
			// default this means 1 closed loop)
			if ( best_count >= (Size)select_best_loop_from ) break;
			// If we can't build this loop due to some major fundamental failure
			if ( result == CriticalFailure ) break;
			// or if we're still unsuccessful (i.e. best_count == 0) and growing the loop isnt an option (because strict_loops is set) then also give up.
			if ( strict_loops_ ) break;
		}


		if ( (best_count > 0) || accept_aborted_loops_ ){
			pose_initial = best_pose;
			pose = best_pose;
		}


		// Print statistics:
		int time_end = time(NULL);
		float time_per_build = float(time_end - time_start) / float(nfailure+1);

		tr.Info   << "Loopstat: "
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
				if ( (Size)loop_numbers[i]  > loops_.size()  ){
					utility_exit_with_message( "Specified loop number is greater than the number of loops in the loop file itself.");
				}
				temp_loops.add_loop( loops_[ loop_numbers[i] ] );
		}

	} else {
		// Choose loops by skiprate
		for ( Loops::const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it ) {
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
		std::random_shuffle( comb_loops.v_begin(), comb_loops.v_end() );
	}

	selected_loops = comb_loops;
} // IndependentLoopMover::select_loops()

std::string
IndependentLoopMover::get_name() const {
	return "IndependentLoopMover";
}

} // namespace loops
} // namespace protocols
