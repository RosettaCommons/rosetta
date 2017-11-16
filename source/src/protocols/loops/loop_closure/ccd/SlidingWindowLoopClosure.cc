// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>

// Package Headers
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/LoopClosure.hh>
#include <protocols/loops/loop_closure/ccd/ShortLoopClosure.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Exceptions.hh>


#include <protocols/evaluation/util.hh>
// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragCache.hh> // for FragStore

#include <core/fragment/SecondaryStructure.hh>

//*only for debug structures
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <basic/options/option.hh> // for quick-test from run:dry_run
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/fast_loops.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

//numeric headers

#include <core/fragment/FragData.hh>
#include <core/fragment/FrameIterator.hh>
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
#include <utility/vector1.hh>

#include <protocols/jobdist/Jobs.hh>
//Auto Headers


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using namespace pose;

static basic::Tracer tr( "protocols.loops.loop_closure.ccd.SlidingWindowLoopClosure" );

std::string const VDW_FRAG_STORE( "closure_loop_vdw" );
std::string const SCORE_FRAG_STORE( "closure_loop_score" );
std::string const RMSD_FRAG_STORE( "closure_loop_rmsd" );
const Real REALLY_BAD_SCORE ( 1000000000.0 );

SlidingWindowLoopClosure::SlidingWindowLoopClosure(
	fragment::FragSetCOP fragset,
	scoring::ScoreFunctionOP scorefxn,
	kinematics::MoveMapCOP movemap
) : scorefxn_( scorefxn ),
	movemap_( movemap ),
	fragset_( fragset ),
	ss_info_( core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure( *fragset ) ) )
{
	set_defaults();
}

SlidingWindowLoopClosure::SlidingWindowLoopClosure() :
	scorefxn_( /* NULL */ ),
	movemap_( /* NULL */ ),
	fragset_( /* NULL */ ),
	ss_info_( /* NULL */ )
{
	set_defaults();
}

SlidingWindowLoopClosure::~SlidingWindowLoopClosure() {}

//@brief sets the movemap
void SlidingWindowLoopClosure::movemap( core::kinematics::MoveMapCOP movemap ) {
	movemap_ = movemap;
}

void
SlidingWindowLoopClosure::set_defaults() {
	using namespace basic::options;

	min_loop_size_ = 6;
	max_loop_size_ = 12;
	min_good_loops_ = 3;
	min_breakout_good_loops_ = 6;
	vdw_delta_ = 0.5;
	score_delta_ = 0.5;
	vdw_score_type_ = scoring::vdw; //not used right now
	bKeepFragments_ = false;
	scored_frag_cycle_ratio_ = 1.0;
	short_frag_cycle_ratio_ = 1.0;
	bQuickTest_ = basic::options::option[ basic::options::OptionKeys::run::dry_run ]();
	bIdealLoopClosing_ = true;
	chainbreak_max_ = 0.2;
	tr.Info << "SlidingWindowLoopClosure::defaults " << std::endl;
}


void SlidingWindowLoopClosure::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		std::string silent_file="_"+prefix;
		if ( get_current_job() && get_current_job()->output_file_name() != "" ) {
			silent_file = get_current_job()->output_file_name()+silent_file;
		} else silent_file = option[ basic::options::OptionKeys::out::file::silent ]+silent_file;

		core::io::silent::SilentFileOptions opts;
		SilentFileData sfd( opts );
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//  ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		pss->fill_struct( pose, get_current_tag() );

		sfd.write_silent_struct( *pss, silent_file, false /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}


//called from apply to figure out which cut is removed, setup min-max loopsizes etc.
Loop
SlidingWindowLoopClosure::determine_loop( Pose const& more_cut, Pose& less_cut ) {
	using namespace kinematics;
	tr.Debug << " SlidingWindowLoopClosure init... " << std::endl;
	FoldTree const& f_more( more_cut.fold_tree() );
	FoldTree const& f_less( less_cut.fold_tree() );
	if ( f_more.num_cutpoint() != f_less.num_cutpoint() + 1 ) {
		utility_exit_with_message(" Failure in SlidingWindowLoopClosure.cc : more_cut should have one more cutpoint then less_cut.  ");
	}

	// which cut is missing ?
	Size cutpoint( 0 );
	for ( Size i = 1; i <= static_cast< Size >( f_more.num_cutpoint() ) && !cutpoint ; i++ ) {
		Size const seqpos( f_more.cutpoint( i ) );
		if ( !f_less.is_cutpoint( seqpos ) ) cutpoint = seqpos;
	}
	//runtime_assert( cutpoint );
	if ( cutpoint == 0 ) { utility_exit_with_message( "Failure in: SlidingWindowLoopClosure.cc: cutpoint seq num is 0 " ); }

	tr.Debug << " remove cutpoint: " << cutpoint << std::endl;
	// find min - max loop
	// figure out where the loop can safely start and stop
	// don't want to have jump_points inside loop because it would create motion downstream of that jump
	//  Size min_loop_begin ( cutpoint + 1 );
	//  Size max_loop_end  ( cutpoint );

	// compute  max_loop_size...
	// extend loop from cutpoint away until either a jump-residue or an unmovable bb-torsion is found
	Size min_loop_begin ( cutpoint + 1 );
	Size max_loop_end  ( cutpoint );
	while (
			min_loop_begin > 1
			&& !f_more.is_jump_point( min_loop_begin - 1 )
			&& movemap().get_bb( min_loop_begin - 1 )
			) --min_loop_begin;
	while (
			max_loop_end < f_more.nres()
			&& !f_more.is_jump_point( max_loop_end + 1 )
			&& movemap().get_bb( max_loop_end + 1 )
			) ++max_loop_end;

	//Size const actual_max_loop_size( max_loop_end_ - min_loop_begin_ + 1 );
	//runtime_assert( min_loop_begin <= max_loop_end );
	if ( min_loop_begin > max_loop_end ) { utility_exit_with_message(" Failure in SlidingWindowLoopClosure.cc : min_loop_begin > max_loop_end "); };

	return Loop( min_loop_begin, max_loop_end, cutpoint );

}


core::scoring::ScoreFunctionOP
SlidingWindowLoopClosure::setup_frag_scorefxn() {
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//   scorefxn.set_weight( vdw, 1.0 );
	//   scorefxn.set_weight( env, 1.0 );
	//   scorefxn.set_weight( pair, 1.0 );
	//   scorefxn.set_weight( linear_chainbreak, 1.0 );
	ScoreFunctionOP score_tmp_;
	tr.Debug << " get fragsample scorefunction... " << std::endl;
	score_tmp_  = ScoreFunctionFactory::create_score_function(
		option[ OptionKeys::fast_loops::fragsample_score ]() /*default "loop_fragsample" );*/
	);
	if ( option[ OptionKeys::fast_loops::fragsample_patch ].user() ) {
		score_tmp_->apply_patch_from_file( option[ OptionKeys::fast_loops::fragsample_patch ] );
	}

	// We will be returning score_tmp_ after we setup the scorefxn_ member variable

	tr.Debug << " get filter scorefunction... " << std::endl;
	if ( option[ OptionKeys::fast_loops::overwrite_filter_scorefxn ].user() ) {
		scorefxn_ = ScoreFunctionFactory::create_score_function( option[ OptionKeys::fast_loops::overwrite_filter_scorefxn ] );
	}

	if ( !scorefxn_ ) {
		scorefxn_ = core::scoring::ScoreFunctionOP( new ScoreFunction );
		scorefxn_->set_weight( linear_chainbreak, 1.0 );
		scorefxn_->set_weight( overlap_chainbreak, 1.0 );
	}
	if ( option[ OptionKeys::fast_loops::patch_filter_scorefxn ].user() ) {
		tr.Debug << " patch filter scorefunction... " << std::endl;
		scorefxn_->apply_patch_from_file( option[ OptionKeys::fast_loops::patch_filter_scorefxn ] );
	}
	if ( scorefxn_->get_weight( linear_chainbreak ) == 0.0 ) {
		scorefxn_->set_weight( linear_chainbreak, 1.0 ); //these have to be switched on...otherwise unclosed loops are invisible
	}
	if ( scorefxn_->get_weight( overlap_chainbreak ) == 0.0 ) {
		scorefxn_->set_weight( overlap_chainbreak, 1.0 );
	}
	tr.Debug << " is there a filter_cst_ evaluator?.. " << std::endl;
	if ( !filter_cst_ && option[ OptionKeys::fast_loops::filter_cst_file ].user() ) {
		filter_cst_ = constraints_additional::ConstraintEvaluatorOP( new constraints_additional::ConstraintEvaluator( "filter_loops", option[ OptionKeys::fast_loops::filter_cst_file ]() ) );
		filter_cst_weight_ = option[ OptionKeys::fast_loops::filter_cst_weight ]();
	}

	return score_tmp_;
}

void SlidingWindowLoopClosure::apply( Pose& more_cut ) {
	// call this only with non-ideal loop closing
	// and with loop_ already set !
	runtime_assert( !bIdealLoopClosing() );
	apply( more_cut, more_cut ); //less_cut loop isn't really used if it is non-ideal loop closing
}

std::string
SlidingWindowLoopClosure::get_name() const {
	return "SlidingWindowLoopClosure";
}


void SlidingWindowLoopClosure::setPoseExtraScore( core::pose::Pose &pose ) {
	core::pose::setPoseExtraScore( pose, "loop_vdw_score", -1 );
	core::pose::setPoseExtraScore( pose, "loop_chain_score", -1 );
	core::pose::setPoseExtraScore( pose, "loop_total_score", -1 );
	core::pose::setPoseExtraScore( pose, "loop_overlap_score", -1 );
	core::pose::setPoseExtraScore( pose, "looprms", -1 );
}

void
SlidingWindowLoopClosure::apply( Pose& more_cut, Pose& less_cut ) {
	if ( closure_fragments() ) { //because we might leave the following block via Exception for some poses
		setPoseExtraScore( more_cut );
	}


	scoring::ScoreFunctionOP frag_scorefxn = setup_frag_scorefxn();
	tr.Debug << "Trying loop-sizes: " << loop_ << std::endl;
	tr.Info << "---------------- LOOP SAMPLING based on this scorefunction: ----------------\n";
	if ( tr.Info.visible() ) frag_scorefxn->show( tr.Info, more_cut );
	tr.Info << std::endl;

	tr.Debug << "Trying loop-sizes: " << loop_ << std::endl;
	tr.Info << "---------------- LOOP SELECTION based on this scorefunction: ----------------\n";
	if ( tr.Info.visible() ) scorefxn_->show( tr.Info, more_cut );
	tr.Info << std::endl;

	loops::remove_cutpoint_variants( more_cut, true );
	loops::add_single_cutpoint_variant( more_cut, loop_ );

	loops::remove_cutpoint_variants( less_cut, true  );
	loops::add_single_cutpoint_variant( less_cut, loop_ );

	tr.Debug << "MOREFOLDTREE: " << more_cut.fold_tree();
	tr.Debug << "LESSFOLDTREE: " << less_cut.fold_tree();
	if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( more_cut, *evaluator_, tr.Debug );
	if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( less_cut, *evaluator_, tr.Debug );


	sample_loops( more_cut, less_cut );
	select_final_loop( more_cut, less_cut );
	loops::remove_cutpoint_variants( more_cut, true  );
	loops::remove_cutpoint_variants( less_cut, true  );
}

void
SlidingWindowLoopClosure::sample_loops( Pose& more_cut, Pose& less_cut ) {
	if ( bQuickTest() ) return; //let's say we found a good loop

	if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ]() ) {
		min_loop_size_ = 6;
		max_loop_size_ = 7;
		min_good_loops_ = 1;
		scored_frag_cycle_ratio_ = 0.05;
		short_frag_cycle_ratio_ = 0.05;
		vdw_delta_ = 200;
		score_delta_ = 200;
		chainbreak_max_ = 20.2;
	}


	best_score_ = REALLY_BAD_SCORE;

	if ( bKeepFragments_ ) {
		closure_fragments_ = core::fragment::OrderedFragSetOP( new fragment::OrderedFragSet );
	}

	Size attempt_count = 0;
	Size const actual_max_loop_size( loop_.size() );
	Size const min_breakout_good_loops( std::max( min_breakout_good_loops_, min_good_loops_ ) );

	tr.Debug << "Trying loop-sizes: "
		<< " "   << std::min( actual_max_loop_size, min_loop_size_ )
		<< " "   << " - "
		<< " "  << std::min( actual_max_loop_size, max_loop_size_ )
		<< " "   << std::min( actual_max_loop_size, max_loop_size_ ) -  std::min( actual_max_loop_size, min_loop_size_ )
		<< " "   << std::endl;

	tr.Debug << "LOOP: " << loop_ << std::endl;

	Size good_loop_count( 0 );

	{ //take initial loop-conformation as closing candidate
		using namespace fragment;
		FrameOP closure_frame( new Frame( loop_.start(), FragDataCOP( FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), loop_.size() ) ) ) ) );
		FrameList closure_frames;
		closure_frame->steal( more_cut );

		CCDLoopClosureMover fast_ccd( loop_, movemap_ );
		fast_ccd.apply( more_cut );
		if ( fast_ccd.success() ) {
			closure_frame->steal( more_cut );
		}

		closure_frames.push_back( closure_frame );

		good_loop_count += process_fragments( closure_frames, more_cut, less_cut );
	} //scope

	scoring::ScoreFunctionOP frag_scorefxn = setup_frag_scorefxn();

	loops::remove_cutpoint_variants( more_cut, true );
	loops::add_single_cutpoint_variant( more_cut, loop_ );

	loops::remove_cutpoint_variants( less_cut, true  );
	loops::add_single_cutpoint_variant( less_cut, loop_ );

	tr.Debug << "MOREFOLDTREE: " << more_cut.fold_tree();
	tr.Debug << "LESSFOLDTREE: " << less_cut.fold_tree();
	if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( more_cut, *evaluator_, tr.Debug );
	if ( evaluator_ && tr.Debug.visible() ) evaluate_pose( less_cut, *evaluator_, tr.Debug );

	// try different loop sizes
	for ( Size loop_size = std::min( actual_max_loop_size, min_loop_size_ );
			loop_size <= std::min( actual_max_loop_size, max_loop_size_ ); ++loop_size ) {
		tr.Debug << "loop-size: " << loop_size << std::endl;
		// try different sliding windows, sorted by loop_fraction in window
		WindowList windows;
		generate_window_list( loop_size, windows );
		windows.sort();
		windows.reverse();

		for ( WindowList::const_iterator it = windows.begin(),
				eit = windows.end(); it != eit; ++it ) {

			attempt_count ++;

			// close this loop
			Loop const current_loop = it->second;
			tr.Debug << "attempt closure on " << current_loop << " # :" << attempt_count << std::endl;
			fragment::FrameList closure_frames;

			tr.Info << "scored fragment sampling... " << current_loop << std::endl;
			LoopClosure scored_closure( fragset_, frag_scorefxn, current_loop, movemap_ );
			scored_closure.set_cycles( scored_frag_cycle_ratio_ );
			scored_closure.ramp_chainbreak();

			//  if ( bIdealLoopClosing() ){
			// note this apply doesn't change the more_cut pose
			if ( scored_closure.apply( more_cut ) ) closure_frames.push_back( scored_closure.closure_fragments() );

			if ( loop_size <= 10 ) {
				tr.Info << "short loop: unscored fragment sampling... " << std::endl;
				ShortLoopClosure short_closure( fragset_, current_loop, movemap_ );
				short_closure.set_cycles( short_frag_cycle_ratio_ );
				short_closure.ramp_chainbreak();
				if ( short_closure.apply( more_cut ) ) closure_frames.push_back( short_closure.closure_fragments() );
			}

			good_loop_count += process_fragments( closure_frames, more_cut, less_cut );
			tr.Info << "process fragments... GoodLoops: " <<  good_loop_count << std::endl;
			if ( good_loop_count >= min_breakout_good_loops ) break;

		} // windows
		if ( good_loop_count >= min_good_loops_ && best_fragment_.is_valid() ) {
			return;
		}
	} // loop_size
	tr.Warning << "no good loop found !" << std::endl;
	throw( loops::EXCN_Loop_not_closed() );
} //apply


void
SlidingWindowLoopClosure::select_final_loop( Pose& more_cut, Pose& less_cut ) {
	using namespace fragment;

	if ( !best_fragment_.is_valid() ) {
		tr.Error << "Cannot apply best fragment because it is not valid or there is none." << std::endl;
		core::pose::setPoseExtraScore( more_cut, "loop_vdw_score", 0 );
		core::pose::setPoseExtraScore( more_cut, "loop_chain_score", 0);
		core::pose::setPoseExtraScore( more_cut, "loop_total_score", 0 );
		core::pose::setPoseExtraScore( more_cut, "loop_overlap_score", 0 );
		core::pose::setPoseExtraScore( more_cut, "looprms", 0 );
		return ;
	}

	best_fragment_.apply( movemap(), less_cut );
	best_fragment_.apply( movemap(), more_cut );

	FragStore< Real > vdw_store(VDW_FRAG_STORE);
	FragStore< Real > score_store(SCORE_FRAG_STORE);
	FragStore< Real > chainbreak_store("chainbreak");
	FragStore< Real > overlap_store("overlap");
	FragStore< Real > rms_store("loop_rms");
	core::pose::setPoseExtraScore( more_cut, "loop_vdw_score", vdw_store.retrieve( best_fragment_ ));
	core::pose::setPoseExtraScore( more_cut, "loop_chain_score", chainbreak_store.retrieve( best_fragment_ ));
	core::pose::setPoseExtraScore( more_cut, "loop_total_score", score_store.retrieve( best_fragment_ ));
	core::pose::setPoseExtraScore( more_cut, "loop_overlap_score", overlap_store.retrieve( best_fragment_ ));
	core::pose::setPoseExtraScore( more_cut, "looprms", rms_store.retrieve( best_fragment_ ));

	if ( closure_fragments() ) {
		std::string frag_file=basic::options::option[ basic::options::OptionKeys::out::file::silent ]()+"_best_frags";
		utility::io::ozstream frag_stream;
		frag_stream.open_append( frag_file );

		for ( FragID_Iterator it = closure_fragments()->begin(),
				eit = closure_fragments()->end(); it != eit; ++it ) {
			Real const size_fraction( 1.0*(*it).frame().length() / loop_.size());
			frag_stream
				<< RJ( 10, vdw_store.retrieve( *it ) )<< " "
				<< RJ( 10, score_store.retrieve( *it ) )<< " "
				<< RJ( 10, chainbreak_store.retrieve( *it) )<< " "
				<< RJ( 10, overlap_store.retrieve( *it ) )<< " "
				<< RJ( 10, rms_store.retrieve( *it ) )<< " "
				<< RJ( 10, size_fraction ) << " "
				<< get_current_tag() << std::endl;
		}
	}
} //select_final_loop

Real
SlidingWindowLoopClosure::filter_score( core::pose::Pose& pose ) {
	Real score ( (*scorefxn_)( pose ) );
	if ( filter_cst_ ) {
		score+=filter_cst_weight_*filter_cst_->apply( pose );
	}
	return score;
}

Size
SlidingWindowLoopClosure::process_fragments(
	fragment::FrameList& frame_in,
	Pose const& more_cut,
	Pose const& loop_pose
) {
	using namespace fragment;
	FragStore< Real > vdw_store(VDW_FRAG_STORE);
	FragStore< Real > score_store(SCORE_FRAG_STORE);
	FragStore< Real > rmsd_store(RMSD_FRAG_STORE);
	FragStore< Real > chainbreak_store("chainbreak");
	FragStore< Real > overlap_store("overlap");
	FragStore< Real > rms_store("loop_rms");

	Size good_loop_count( 0 );

	Loops loops;
	loops.add_loop( loop_ ); //for looprms

	Pose orig_pose;
	Pose work_pose;

	if ( bIdealLoopClosing() ) {
		orig_pose = more_cut;
		work_pose = loop_pose;
	} else {
		orig_pose = more_cut;
		work_pose = more_cut;
	}

	Size ct( 0 );

	Real const start_score ( filter_score( orig_pose ) );
	Real const vdw_start_score ( orig_pose.energies().total_energies()[ vdw_score_type_ ] );
	for ( FragID_Iterator it = frame_in.begin(),
			eit = frame_in.end(); it != eit; ++it ) {
		if ( !it->is_valid() ) continue;
		if ( it->apply( movemap(), work_pose ) ) {
			ct++;

			//check rmsd
			if ( bIdealLoopClosing() && !basic::options::option[ basic::options::OptionKeys::run::test_cycles ]() ) {
				it->apply( movemap(), orig_pose );
				Real const check_rmsd( scoring::CA_rmsd( orig_pose, work_pose ) ); //, loop_.start(), loop_.stop() ) );
				// Real sssscore = (*scorefxn_)( orig_pose );
				// Real const lin_cb( orig_pose.energies().total_energies()[ scoring::linear_chainbreak ]);
				// Real const ovp_cb( orig_pose.energies().total_energies()[ scoring::overlap_chainbreak ]);
				//   tr.Debug << "CHECK_LOOP: " << check_rmsd << " " << lin_cb << " " << ovp_cb << " " << sssscore << std::endl;
				if ( check_rmsd > 0.1 ) {
					tr.Debug << "skip loop: bad check_rms in loop building RMSD: "<< check_rmsd << std::endl;
					//throw them out...
					continue;
				}
			}

			Real const score( filter_score( work_pose ) );
			Real const vdw_score( work_pose.energies().total_energies()[ vdw_score_type_ ]);
			Real const chainbreak_scores( work_pose.energies().total_energies()[ scoring::linear_chainbreak ]
				+ work_pose.energies().total_energies()[ scoring::chainbreak ]
				+ work_pose.energies().total_energies()[ scoring::overlap_chainbreak ]
			);
			Real const loop_rmsd( get_native_pose() ? loops::loop_rmsd( *get_native_pose(), work_pose, loops, true /*CA*/ ) : -1.0 );
			//   scorefxn_->show( tr.Debug, work_pose ); tr.Trace<< std::endl;
			bool const good_loop( ( vdw_score <= vdw_start_score + vdw_delta_ ) && chainbreak_scores <= chainbreak_max_ );
			//   bool const good_loop( score < start_score + score_delta_ );


			if ( good_loop ) { //originally this was vdw_score < vdw_start_score... is this better?!
				tr.Info << "good loop "<< ct << " "  << score << " " << start_score << " vdw: " << vdw_score << " " << vdw_start_score << std::endl;
				++good_loop_count;
				vdw_store.store( *it, vdw_score-vdw_start_score);
				score_store.store( *it, score );
				chainbreak_store.store( *it, chainbreak_scores );
				overlap_store.store( *it, work_pose.energies().total_energies()[ scoring::overlap_chainbreak ] );
				rms_store.store( *it, loop_rmsd );
				//store good fragment if we are set up to do so
				if ( closure_fragments_ ) {
					// keep these data with the fragments
					closure_fragments_->add( *it );
				}

				if ( best_score_ > score ) { //CHANGED 3/12/09 before we would also take fragments as "best" if they are not "good"... seems stupid...     tr.Debug << "improved loop score: " << score << " prev best: " << best_score_  << std::endl;
					if ( !good_loop ) { //
						tr.Debug << "...but on bad loop " << score << " " << start_score << " vdw: " << vdw_score << " " << vdw_start_score
							<< "chainbreak: " << chainbreak_scores << std::endl;
					}
					best_score_ = score;
					best_fragment_ = *it;
				}
				tr.Debug << ct << " LoopScore: " << score << " Best: " << best_score_ << std::endl;
			} else {
				tr.Debug << ct << " discard loop-fragment ... score: " << score << " chainbreak_scores: " << chainbreak_scores << std::endl; // vdw < vdw_start_score_
			}
		} // an applicable fragment
	} // FragID iteration
	return good_loop_count;
}

void
SlidingWindowLoopClosure::generate_window_list( Size loop_size, WindowList& window_list ) const {
	runtime_assert( ss_info_ != 0 );
	for ( Size ii= std::max( 0, - (int) loop_.cut() + (int) loop_size )  ; ii<= loop_size; ++ii ) {
		runtime_assert( (int) loop_.cut() - (int) loop_size + ii + 1  > 0 );
		Size const loop_begin( loop_.cut() - loop_size + ii + 1);
		Size const loop_end  ( loop_.cut() + ii );
		tr.Debug << "add-window: " << loop_begin << "-" << loop_end << std::endl;
		runtime_assert( loop_begin <= loop_.cut()+1);
		runtime_assert( loop_end >= loop_.cut() );
		runtime_assert( loop_end == loop_begin + loop_size - 1 );
		if ( loop_begin < loop_.start() || loop_end > loop_.stop() ) continue;
		Real f(0);
		for ( Size i = loop_begin; i <= loop_end; ++i ) {
			f += ss_info_->loop_fraction(i);
		}
		window_list.push_back( std::make_pair( f/(1.0*loop_size)-loop_size*1.0, Loop( loop_begin, loop_end, loop_.cut() ) ) );
	}
}

void
SlidingWindowLoopClosure::fragments( core::fragment::FragSetCOP frags ) {
	ss_info_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure( *frags ) );
	fragset_ = frags;
	debug_assert( loop_.start() == 0 || fragset_->max_pos() > (loop_.start() + std::min( loop_.size(), max_loop_size_ )) );
}

void
SlidingWindowLoopClosure::set_loop( Loop const& loop_in ) {
	loop_ = loop_in;
	if ( fragset_ &&
			fragset_->max_pos() <= loop_.start() + std::min( loop_.size(), max_loop_size_ ) ) {
		std::ostringstream ss;
		ss << this->get_name() << " was given a loop to close (" << loop_
			<< ") for which it does not have any fragments (residues: " << fragset_->min_pos() << "-"
			<< fragset_->max_pos() << ").";
		throw utility::excn::EXCN_Msg_Exception( ss.str() );
	}
}


} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
