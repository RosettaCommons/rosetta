// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/loops/loop_closure/ccd/FASelectSlidingWindowLoopClosure.hh>

// Package Headers
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Exceptions.hh>
#include <protocols/loops/util.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragCache.hh> // for FragStore


//#include <protocols/jumping/SecondaryStructure.hh>

#include <basic/options/option.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/fast_loops.OptionKeys.gen.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
//numeric headers

#include <basic/options/option_macros.hh>

#include <core/fragment/FragData.hh>
#include <core/fragment/FrameIterator.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/relax/FastRelax.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


OPT_1GRP_KEY( Real, fast_loops, rmsd_dump )

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using namespace pose;

static thread_local basic::Tracer tr( "protocols.loops.loop_closure.ccd.FASelectSlidingWindowLoopClosure" );

const Real REALLY_BAD_SCORE ( 1000000000.0 );

FASelectSlidingWindowLoopClosure::FASelectSlidingWindowLoopClosure(
  fragment::FragSetCOP fragset,
  scoring::ScoreFunctionOP scorefxn,
  kinematics::MoveMapCOP movemap
) : WidthFirstSlidingWindowLoopClosure( fragset, scorefxn, movemap )
{
	set_defaults();
}

FASelectSlidingWindowLoopClosure::FASelectSlidingWindowLoopClosure() :
	WidthFirstSlidingWindowLoopClosure()
{
	set_defaults();
}

FASelectSlidingWindowLoopClosure::~FASelectSlidingWindowLoopClosure() {}

std::string
FASelectSlidingWindowLoopClosure::get_name() const {
	return "FASelectSlidingWindowLoopClosure";
}

void
FASelectSlidingWindowLoopClosure::set_defaults() {
	Parent::set_defaults();
	keep_fragments(); //we need this... not only for debug
}

void FASelectSlidingWindowLoopClosure::register_options() {
  using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( fast_loops::rmsd_dump, "dump all pdbs for loops that are below in rmsd", 2.0 );
}


void
FASelectSlidingWindowLoopClosure::select_final_loop( Pose& more_cut, Pose& less_cut ) {
	using namespace fragment;
	runtime_assert( closure_fragments() != 0 );
	// in case that we jump out... add extra scrores
	core::pose::setPoseExtraScore( more_cut, "post_relax_looprms", -1 );
  core::pose::setPoseExtraScore( more_cut, "fa_score", -1 );

	Loops loops;
	loops.add_loop( loop_ ); //for looprms


		//sort list of fragments by score
	typedef std::list< std::pair< core::Real, FragID > > FragList;
	FragList scored_frags;

 	FragStore< Real > score_store(SCORE_FRAG_STORE);
	FragStore< Real > rms_store("loop_rms");
	FragStore< Real > fascore_store("fa_score");
	FragStore< Real > post_relax_rms_store("post_relax_looprms");
	FragStore< Real > vdw_store(VDW_FRAG_STORE);
	FragStore< Real > chainbreak_store("chainbreak");
	FragStore< Real > overlap_store("overlap");


	for ( FragID_Iterator it = closure_fragments()->begin(),
					eit = closure_fragments()->end(); it != eit; ++it ) {
		scored_frags.push_back( std::make_pair( score_store.retrieve( *it ), *it ) );
	}
	scored_frags.sort();
	scored_frags.reverse();
	core::Size const Ntot( scored_frags.size() );
	core::Size const Ntest( Ntot ); //0.1*Ntot+0.5 );
	tr.Debug << "try " << Ntest << " of " <<  Ntot << " fragments in full-atom mode" << std::endl;

	best_score_ = REALLY_BAD_SCORE;
	FragID best_fragment_;

	///* score 10% frags with full-atom
	std::string frag_file=basic::options::option[ basic::options::OptionKeys::out::file::silent ]()+"_best_frags";

	FragList::iterator frag_list_it = scored_frags.begin();
	for ( Size ct = 1; ct <= Ntest; ct ++, ++frag_list_it ) {
		FragID const& frag( frag_list_it->second );

		Pose centroid_pose( more_cut );
		//fragments are smaller than the full loop region ...
		// apply to more_cut and steal full loop region to transfer to fa_pose
		frag.apply( movemap(), centroid_pose );

		Pose fa_pose( *fa_pose_ );
		fa_pose.fold_tree( more_cut.fold_tree() );
		set_extended_torsions_and_idealize_loops( fa_pose, loops );
		add_single_cutpoint_variant( fa_pose, loop_ );

		if ( tr.Trace.visible() )	fa_pose.dump_pdb("set_fold_tree_fa_pose.pdb");


		FrameOP fa_frame_( new Frame( loop_.start(), FragDataCOP( FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), loop_.size() ) ) ) ) );
		fa_frame_->steal( centroid_pose );
		fa_frame_->apply( 1, fa_pose );

		if ( tr.Trace.visible() )		fa_pose.dump_pdb("fa_pose.pdb");
		if ( tr.Trace.visible() )		centroid_pose.dump_pdb( "centroid.pdb" );

		if ( rms_store.retrieve( frag ) < basic::options::option[ basic::options::OptionKeys::fast_loops::rmsd_dump ]() ) {
			output_debug_structure( centroid_pose, "_centroid");
			output_debug_structure( fa_pose, "_fa_pre_relax");
		}

		// full-atom relax
		Real fa_score( fascore( fa_pose ) );

		if ( rms_store.retrieve( frag ) < basic::options::option[ basic::options::OptionKeys::fast_loops::rmsd_dump ]() ) {
			output_debug_structure( fa_pose, "_fa_post_relax");
		}

		Real post_relax_looprms( get_native_pose() ? loops::loop_rmsd( *get_native_pose(), fa_pose, loops, true /*CA*/ ) : -1.0 );
		if ( tr.Trace.visible() )		fa_pose.dump_pdb("fa_relaxed_pose.pdb");

		if ( fa_score < best_score_ ) {
			best_score_ = fa_score;
			best_fragment_ = frag;
		}

		fascore_store.store( frag, fa_score );
	  post_relax_rms_store.store( frag, post_relax_looprms );

		utility::io::ozstream frag_stream;
		frag_stream.open_append( frag_file );

		Real const size_fraction( 1.0*frag.frame().length() / loop_.size());
		frag_stream << RJ( 10, vdw_store.retrieve( frag ) )<< " "
								<< RJ( 10, score_store.retrieve( frag ) )<< " "
								<< RJ( 10, chainbreak_store.retrieve( frag) )<< " "
								<< RJ( 10, overlap_store.retrieve( frag ) )<< " "
								<< RJ( 10, rms_store.retrieve( frag ) )<< " "
								<< RJ( 10, fascore_store.retrieve( frag ) )<< " "
								<< RJ( 10, post_relax_rms_store.retrieve( frag ) )<< " "
								<< RJ( 10, size_fraction ) << " "
								<< get_current_tag() << std::endl;
	}


	if ( ! best_fragment_.frame().is_valid() ) {
		throw( loops::EXCN_Loop_not_closed() );
	}
	//* apply winner and generate output, --- essentially copy from Parent class
	best_fragment_.apply( movemap(), less_cut );
  best_fragment_.apply( movemap(), more_cut );

  core::pose::setPoseExtraScore( more_cut, "loop_vdw_score", vdw_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "loop_chain_score", chainbreak_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "loop_total_score", score_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "loop_overlap_score", overlap_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "looprms", rms_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "post_relax_looprms", post_relax_rms_store.retrieve( best_fragment_ ));
  core::pose::setPoseExtraScore( more_cut, "fa_score", fascore_store.retrieve( best_fragment_ ));


}
 //select_final_loop

core::Real FASelectSlidingWindowLoopClosure::fascore( Pose& pose ) const {
	using namespace basic::options;
	Loops loops;
	loops.add_loop( loop_ ); //for looprms

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	scorefxn->show( std::cout, pose );
	std::cout << std::endl;
	// default move map
	(*scorefxn)( pose ); //score to update tenAneighbour graph
	bool const fix_natsc = basic::options::option[ basic::options::OptionKeys::loops::fix_natsc ];
	kinematics::MoveMapOP mm_all_loops( new kinematics::MoveMap ); // DJM tmp
	loops_set_move_map( pose, loops, fix_natsc, *mm_all_loops); // DJM tmp

	//	scorefxn->set_weight( scoring::chainbreak, 10 );
	relax::FastRelax fast_relax( scorefxn, option[ OptionKeys::fast_loops::fast_relax_sequence_file ]() );
	fast_relax.set_movemap( mm_all_loops );
	fast_relax.apply( pose );
	return (*scorefxn)( pose );
}

void FASelectSlidingWindowLoopClosure::set_fullatom_pose( core::pose::Pose& fa_pose ) {
	fa_pose_ = core::pose::PoseOP( new core::pose::Pose( fa_pose ) );
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
