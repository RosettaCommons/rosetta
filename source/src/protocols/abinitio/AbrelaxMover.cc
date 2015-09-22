// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/abinitio/AbrelaxMover.cc
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

// Unit Headers
#include <protocols/abinitio/AbrelaxMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/abinitio/FragmentSampler.hh>
#include <protocols/abinitio/ConstraintFragmentSampler.hh>
//#include <protocols/abinitio/TMHTopologySampler.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jumping/util.hh>
#include <protocols/loops/Exceptions.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/util.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/simple_moves/RepulsiveOnlyMover.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/abrelax.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <numeric/random/random.hh>

#include <utility/vector1.hh>

#include <protocols/moves/PyMolMover.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.general_abinitio", basic::t_info );


namespace protocols {
namespace abinitio {

using namespace core;

AbrelaxMover::AbrelaxMover() :
	topology_broker_( /* NULL */ ),
	sampling_protocol_( /* NULL */ ),
	loop_closure_protocol_( /* NULL */ ),
	relax_protocol_( /* NULL */ ),
	post_loop_closure_protocol_( /* NULL */ ),
	b_return_unrelaxed_fullatom_( false )
{
	basic::mem_tr << "AbrelaxMover CStor start" << std::endl;
	set_defaults();
	basic::mem_tr << "AbrelaxMover CStor end" << std::endl;
}

AbrelaxMover::~AbrelaxMover() {}

void AbrelaxMover::clear() {
	topology_broker_ = NULL;
	sampling_protocol_ = NULL;
	loop_closure_protocol_ = NULL;
	relax_protocol_ = NULL;
	post_loop_closure_protocol_ = NULL;
}

void AbrelaxMover::set_defaults() {
	using namespace basic::options;
	using protocols::idealize::IdealizeMover;
	using protocols::idealize::IdealizeMoverOP;
	using protocols::loops::loop_closure::ccd::SlidingWindowLoopClosure;
	using protocols::loops::loop_closure::ccd::SlidingWindowLoopClosureOP;
	using protocols::loops::loop_closure::ccd::WidthFirstSlidingWindowLoopClosure;

	b_return_unrelaxed_fullatom_ = false;

	topology_broker_ = topology_broker::TopologyBrokerOP( new topology_broker::TopologyBroker() );
	topology_broker::add_cmdline_claims( *topology_broker_ );


	// FragmentSamplerOP sampler;
	// use the TMHTopologySampler or default ConstraintFragmentSampler
	// if(option[OptionKeys::abinitio::TMH_topology].user())
	// {
	//  tr << "setting TMYTopologySampler" << std::endl;
	//  FragmentSamplerOP sampler = new TMHTopologySampler(topology_broker_);
	//  tr << "sampler:  " << sampler->get_name() << std::endl;
	//  sampling_protocol( sampler );
	// }
	// else{
	//tr << "setting ConstraintfragmentSampler" << std::endl;
	FragmentSamplerOP sampler( new ConstraintFragmentSampler( topology_broker_ ) );
	//tr << "sampler:  " << sampler->get_name() << std::endl;
	sampling_protocol( sampler );
	// }


	//  Idealize the structure before relax
	bool bIdeal( true );
	if ( option[ OptionKeys::loops::idealize_before_loop_close ].user() ) {
		IdealizeMoverOP idealizer( new IdealizeMover );
		idealizer->fast( false );
		pre_loop_closure_protocol( idealizer );
	}

	// loop closing
	if ( option[ OptionKeys::abinitio::close_loops ]() ) {
		SlidingWindowLoopClosureOP closure_method( new SlidingWindowLoopClosure );

		if ( option[ OptionKeys::loops::alternative_closure_protocol ]() ) {
			closure_method = SlidingWindowLoopClosureOP( new WidthFirstSlidingWindowLoopClosure );
		}

		bIdeal = !option[ OptionKeys::loops::non_ideal_loop_closing ]();

		// set options here if you like
		// closure_protocol-> ... and write the setters/accessors, too, if you have to
		closure_method->scored_frag_cycle_ratio( option[ OptionKeys::loops::scored_frag_cycles ]() );
		closure_method->short_frag_cycle_ratio( option[ OptionKeys::loops::short_frag_cycles ]() );
		closure_method->set_bIdealLoopClosing( bIdeal );
		closure_method->set_chainbreak_max( option[ OptionKeys::loops::chainbreak_max_accept ]() );

		//add closure protocol to abrelax
		closure_protocol( closure_method );
	}

	//  Idealize the structure after relax
	if ( option[ OptionKeys::loops::idealize_after_loop_close ].user() ) {
		IdealizeMoverOP idealizer( new IdealizeMover );
		idealizer->fast( false );
		post_loop_closure_protocol( idealizer );
		bIdeal = true;
	}

	//yipee a bad hack using global variables
	if ( !bIdeal ) option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary");

	//add relax protocol to abrelax
	//cst_fa_weight and cst_weight are properly interpreted by this utility function
	relax_protocol( relax::generate_relax_from_cmd( true /*null if no relax flag*/ ) );

	if ( relax_protocol() && option[ OptionKeys::run::test_cycles ] ) {
		//in test_cycles we don't test relax...
		//special relax integration test takes care of that...
		// so make it fast here... dry_run means just an energy evaluation
		relax_protocol()->dry_run( true );
	}
}

FragmentSamplerOP AbrelaxMover::sampling_protocol() {
	return sampling_protocol_;
}

relax::RelaxProtocolBaseCOP AbrelaxMover::relax_protocol() const {
	return relax_protocol_;
}

relax::RelaxProtocolBaseOP AbrelaxMover::relax_protocol() {
	return relax_protocol_;
}

loops::loop_closure::ccd::SlidingWindowLoopClosureOP AbrelaxMover::closure_protocol() {
	return loop_closure_protocol_;
}

void AbrelaxMover::sampling_protocol( FragmentSamplerOP set) {
	sampling_protocol_ = set;
}

void AbrelaxMover::relax_protocol( relax::RelaxProtocolBaseOP set ) {
	relax_protocol_ = set;
}

void AbrelaxMover::closure_protocol( loops::loop_closure::ccd::SlidingWindowLoopClosureOP set ) {
	loop_closure_protocol_ = set;
}

// idealize after loop-closing
void AbrelaxMover::post_loop_closure_protocol( moves::MoverOP move ) {
	post_loop_closure_protocol_ = move;
}

// idealize before loop-closing
void AbrelaxMover::pre_loop_closure_protocol( moves::MoverOP move ) {
	pre_loop_closure_protocol_ = move;
}

topology_broker::TopologyBrokerOP AbrelaxMover::topology_broker() {
	return topology_broker_;
}

//@brief basic apply for the generalized protocol:
void AbrelaxMover::apply( pose::Pose &pose ) {
	// Can we add the PyMOL mover here?
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[OptionKeys::run::show_simulation_in_pymol].user()
				&& option[OptionKeys::run::show_simulation_in_pymol].value() > 0.0 ) {
			protocols::moves::AddPyMolObserver(pose,
				option[OptionKeys::run::keep_pymol_simulation_history](),
				option[OptionKeys::run::show_simulation_in_pymol].value());
		}
	}

	using namespace basic::options;

	runtime_assert( sampling_protocol() != 0 );
	runtime_assert( topology_broker() != 0 );

	tr.Info << "AbrelaxMover: " << get_current_tag() << std::endl;
	basic::mem_tr << "AbrelaxMover::apply" << std::endl;
	basic::show_time( tr,  "AbrelaxMover: start..."+jd2::current_batch()+" "+jd2::current_output_name() );

	// kidnap sampling_protocols's checkpointer - this could ultimately be a singleton i guess
	checkpoint::CheckPointer &checkpoints = sampling_protocol()->get_checkpoints();

	// need to save numeric::random::rg() states such that choices for constraints and fold-tree are the same.
	if ( ! checkpoints.recover_checkpoint( pose, get_current_tag(), "rg_state") ) {
		checkpoints.checkpoint( pose, get_current_tag(), "rg_state");
	}

	checkpoints.debug( get_current_tag(), "rg_state", numeric::random::rg().uniform() );
	if ( sampling_protocol() ) sampling_protocol()->set_current_tag( get_current_tag() );
	if ( closure_protocol() )  closure_protocol()->set_current_tag( get_current_tag() );
	if ( relax_protocol() )    relax_protocol()->set_current_tag( get_current_tag() );

	//   setup -- e.g., sequence --> pose
	//
	topology_broker()->apply( pose ); //creates pose and a state in the topology_broker needed for the whole run

	// apply a mover which calculates only repulsive energy on designate residues
	{
		protocols::simple_moves::RepulsiveOnlyMover replonly;
		replonly.set_mutate_to_glycine( false );
		replonly.apply( pose );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );  //add viewer

#ifdef BOINC_GRAPHICS
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

	//   abinitio
	//
	sampling_protocol()->topology_broker( topology_broker() );

#ifdef BOINC_GRAPHICS
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

	sampling_protocol()->apply( pose );
	if ( sampling_protocol()->get_last_move_status() == protocols::moves::FAIL_RETRY ) {
		this->set_last_move_status(protocols::moves::FAIL_RETRY);
		return;
	}

	//make sure all chainbreak variants are activated:
	if ( !topology_broker()->check_chainbreak_variants( pose ) ) {
		tr.Warning << "[WARNING] some chainbreaks in " << jd2::current_output_name()
			<< " were not penalized at end of " << sampling_protocol()->type()
			<< std::endl;
		topology_broker()->add_chainbreak_variants( pose );
	}

	//   filters
	//
	bool loop_success = true;

#ifdef BOINC_GRAPHICS
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

	scoring::ScoreFunctionCOP last_scorefxn_cop(NULL); // holding handle for clone
	scoring::ScoreFunction const * last_scorefxn ( & sampling_protocol()->current_scorefxn() );

	// Make sure score columns always the same.
	relax::ClassicRelax().setPoseExtraScore( pose );
	protocols::loops::loop_closure::ccd::SlidingWindowLoopClosure::setPoseExtraScore( pose );

	if ( pre_loop_closure_protocol_ ) {
		tr.Info << "abrelax_stage: pre_loopclosing (i.e., idealize) for " << get_current_tag() << std::endl;
		kinematics::FoldTree fold_tree( pose.fold_tree() );
		pose.fold_tree( topology_broker()->final_fold_tree() );
		pre_loop_closure_protocol_->apply( pose );
		pose.fold_tree( fold_tree );
		sampling_protocol()->current_scorefxn()( pose );
		jd2::output_intermediate_pose( pose, "loops_closed_preprocessed" );
	}

	//   loop closing
	//
	if ( closure_protocol() ) {
		// Did we already close the loops successfully ?
		if ( checkpoints.recover_checkpoint( pose, get_current_tag(), "loops_S", false /*fullatom*/, true /*foldtree*/) )  {
			checkpoints.debug( get_current_tag(), "close_loops", (*last_scorefxn)(pose), (core::Real) true );
			loop_success = true;
		} else if ( checkpoints.recover_checkpoint( pose, get_current_tag(), "loops_C", false /*fullatom*/, true /*foldtree*/) )  {
			// No ? Have we already tried but failed ?
			checkpoints.debug( get_current_tag(), "close_loops", (*last_scorefxn)(pose), (core::Real) false );
			loop_success = false;
		} else {
			// No ? Then evidently we havn't even tried yet
			tr << "AbrelaxMover: start loops" << std::endl;
			kinematics::MoveMapOP movemap( new kinematics::MoveMap );
			closure_protocol()->scorefxn( sampling_protocol()->current_scorefxn().clone() );
			closure_protocol()->fragments( topology_broker()->loop_frags( *movemap ) ); //get frags and movemap from broker
			closure_protocol()->movemap( movemap ); //this is the movemap sanctioned by the broker ... see line above

			try {
				jumping::close_chainbreaks(closure_protocol(),
					pose,
					sampling_protocol()->get_checkpoints(),
					get_current_tag(),
					topology_broker()->final_fold_tree());
			} catch ( loops::EXCN_Loop_not_closed& excn ) {
				set_current_tag( "C_"+get_current_tag().substr(std::min(2,(int)get_current_tag().size())) );
				set_last_move_status( moves::FAIL_DO_NOT_RETRY );
				loop_success = false;
				if ( option[ OptionKeys::abrelax::fail_unclosed ]() ) return;
			}
			sampling_protocol()->current_scorefxn()( pose );
			jd2::output_intermediate_pose( pose, "loops_closed" );

			topology_broker()->apply_filter( pose, LOOP_CLOSURE, 1 );

			if ( post_loop_closure_protocol_ && loop_success ) {
				tr.Info << "abrelax_stage: post_loopclosing (i.e., idealize) for " << get_current_tag() << std::endl;
				post_loop_closure_protocol_->apply( pose );
				sampling_protocol()->current_scorefxn()( pose );
				jd2::output_intermediate_pose( pose, "loops_closed_postprocessed" );
			}

			// to know this we'd have to catch the Exception EXCN_Loop_not_closed
			if ( loop_success ) {
				checkpoints.checkpoint( pose, get_current_tag(), "loops_S", true /*foldtree*/ );
			} else {
				checkpoints.checkpoint( pose, get_current_tag(), "loops_C", true /*foldtree*/ );
			}

			checkpoints.debug( get_current_tag(), "loops", (*last_scorefxn)(pose), (core::Real) loop_success );
		}
	}

	//   Fullatom switch
	//
	basic::mem_tr << "AbrelaxMover::apply fullatom switch" << std::endl;
	if ( relax_protocol() || b_return_unrelaxed_fullatom_ ) {
		tr << "AbrelaxMover: switch to fullatom" << std::endl;
		jd2::get_current_job()->add_string_real_pair( "prefa_centroid_score",  ((sampling_protocol()->current_scorefxn())( pose ) ) );
		tr.Info << "prefa_centroid_score:\n ";
		(sampling_protocol()->current_scorefxn()).show( tr.Info, pose );
		tr.Info << std::endl;
		core::scoring::ScoreFunctionOP clean_score3( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ) );
		clean_score3->set_weight( scoring::linear_chainbreak, 1.33 );
		clean_score3->set_weight( scoring::overlap_chainbreak, 1.0 );
		//  clean_score3->set_weight( scoring::chainbreak, 1.0 );
		jd2::get_current_job()->add_string_real_pair( "prefa_clean_score3", ((*clean_score3)(pose)) );
		tr.Info << "prefa_clean_score3:\n ";
		clean_score3->show( tr.Info, pose );
		tr.Info << std::endl;
		topology_broker()->switch_to_fullatom( pose );

		// we're now in fullatom mode - so upate the score function.
		if ( ( (option[ OptionKeys::score::weights ]() == "score0") ||
				(option[ OptionKeys::score::weights ]() == "score2") ||
				(option[ OptionKeys::score::weights ]() == "score3") ||
				(option[ OptionKeys::score::weights ]() == "score5") ) ) {
			utility_exit_with_message("Cannot proceed - you chose a centroid score function for fullatom mode");
		}

		protocols::relax::RelaxProtocolBaseCOP relax_prot = relax_protocol();
		last_scorefxn_cop = relax_prot->get_scorefxn();
		last_scorefxn = last_scorefxn_cop.get();

	} //switched to fullatom


	{// apply a mover which calculates only repulsive energy on designate residues
		protocols::simple_moves::RepulsiveOnlyMover replonly;
		replonly.set_mutate_to_glycine( false );
		replonly.apply( pose );
	}

	if ( option[ basic::options::OptionKeys::abinitio::close_loops_by_idealizing ]() ) {
		close_with_idealization( pose );
	}

	//   Relax
	//
	basic::mem_tr << "AbrelaxMover::apply relax" << std::endl;
	if ( ( loop_success || option[ OptionKeys::abinitio::relax_failures ]() ) && relax_protocol() ) {
		tr << "AbrelaxMover: relax " << std::endl;
		if ( !checkpoints.recover_checkpoint( pose, get_current_tag(), "relax", true, true) ) {
			topology_broker()->adjust_relax_movemap( *(relax_protocol()->get_movemap() ) );
			relax_protocol()->apply( pose );
			checkpoints.checkpoint( pose, get_current_tag(), "relax", true ); //since relax_protocol throws away its checkpoints right here
		}
		checkpoints.debug( get_current_tag(), "relax", (*last_scorefxn)( pose ) );
	}

	topology_broker()->apply_filter( pose, RELAX, 1 );

	//   Final clean relax
	//
	if ( ( loop_success || option[ OptionKeys::abinitio::relax_failures ]() ) && relax_protocol() && option[ OptionKeys::abinitio::final_clean_relax ]() ) {
		tr << "AbrelaxMover: final relax " << std::endl;
		pose.constraint_set( NULL );
		if ( !checkpoints.recover_checkpoint( pose, get_current_tag(), "finalrelax", true, true) ) {
			relax_protocol()->apply( pose );
			checkpoints.checkpoint( pose, get_current_tag(), "finalrelax", true ); //since relax_protocol throws away its checkpoints right here
		}
		checkpoints.debug( get_current_tag(), "finalrelax", (*last_scorefxn)( pose ) );
	}

	if ( option[ OptionKeys::abinitio::clear_pose_cache ]() ) {
		tr.Debug << "\n******************************************************** \n"
			<< "              CLEAR POSE CACHE                             \n"
			<< "*************************************************************"
			<< std::endl;
		pose.data().clear();
	}

	if ( sampling_protocol() ) sampling_protocol()->get_checkpoints().clear_checkpoints();
	if ( !b_return_unrelaxed_fullatom_ ) ( *last_scorefxn)( pose );

	basic::mem_tr << "AbrelaxMover::apply end" << std::endl;
	basic::show_time( tr,  "AbrelaxMover: finished ..."+jd2::current_batch()+" "+jd2::current_output_name() );
}

std::string AbrelaxMover::get_name() const {
	return "AbrelaxMover";
}

void AbrelaxMover::close_with_idealization( pose::Pose &pose) {
	using namespace basic::options;
	kinematics::FoldTree fold_tree( pose.fold_tree() );

	// record cutpoints
	protocols::loops::LoopsOP cloops( new protocols::loops::Loops() );
	for ( Size ncut = 1; ncut <= (Size) pose.fold_tree().num_cutpoint(); ncut++ ) {
		Size cutpoint = pose.fold_tree().cutpoint( ncut );
		Size margin = option[ basic::options::OptionKeys::abinitio::optimize_cutpoints_margin ]();
		protocols::loops::Loop newloop (std::max( (int) 1, int(cutpoint - margin) ),
			std::min( (int) pose.total_residue(), int(cutpoint + margin) ),
			0);
		if ( cloops->size() >= 2 ) {
			if ( newloop.start() <= ( *cloops )[cloops->size()-1].stop() ) newloop.set_start( ( *cloops )[cloops->size()-1].stop() +2 );
		}
		newloop.choose_cutpoint( pose );
		cloops->add_loop( newloop );
	}
	cloops->auto_choose_cutpoints( pose );

	// remove cuts
	pose.fold_tree( topology_broker()->final_fold_tree() );

	//idealize
	protocols::idealize::IdealizeMover idealizer;
	idealizer.fast( false );
	idealizer.apply( pose );

	// Relax, with none of the cutpoints
	relax_protocol()->apply( pose );

	// now optimize around the old cut points.
	if ( option[ basic::options::OptionKeys::abinitio::optimize_cutpoints_using_kic ]() ) {
		protocols::loops::fold_tree_from_loops( pose, *cloops, fold_tree , true /* include terminal cutpoints */);
		pose.fold_tree( fold_tree );

		protocols::relax::RelaxProtocolBaseCOP relax_prot = relax_protocol();
		core::scoring::ScoreFunctionOP refine_scorefxn = relax_prot->get_scorefxn()->clone();
		protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine_kic( cloops, refine_scorefxn );
		refine_kic.apply( pose );

		// Return fold tree to normal state
		pose.fold_tree( topology_broker()->final_fold_tree() );
	}
}

}
}
