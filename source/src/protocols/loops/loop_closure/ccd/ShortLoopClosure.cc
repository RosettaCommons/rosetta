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
#include <protocols/loops/loop_closure/ccd/ShortLoopClosure.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

#include <core/chemical/ChemicalManager.hh>

#include <protocols/simple_moves/FragmentMover.fwd.hh>


#include <protocols/moves/MonteCarlo.hh>
// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
//numeric headers

//// C++ headers
#include <string>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

//Auto Headers


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using namespace pose;

static THREAD_LOCAL basic::Tracer tr( "protocols.loops.loop_closure.ccd.ShortLoopClosure" );


//c'stor
ShortLoopClosure::ShortLoopClosure(
	fragment::FragSetCOP fragset,
	Loop loop_def,
	kinematics::MoveMapCOP
) : LoopClosure(),
	orig_loop_( loop_def )
{
	scoring::ScoreFunctionOP scorefxn( new scoring::ScoreFunction );
	scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
	scorefxn->set_weight( scoring::overlap_chainbreak, 3.0 );
	set_scorefxn( scorefxn );

	set_temperature( 0.0 );

	//  this is now done in base-class
	kinematics::MoveMapOP movemap( new kinematics::MoveMap );
	//   //*movemap = *movemap_in; no copying since seqpos don't make sense
	movemap->set_bb( true ); //alternatively go through "movemap_in and translate settings for each bb
	movemap->set_bb( 1, false );
	movemap->set_bb( loop_def.size()+2, false); //two padding residues on each side of loop
	set_movemap( movemap );

	runtime_assert( loop_def.start() >= 2 );
	int const loop_offset(loop_def.start() - 2); //offset
	Loop short_loop( 2,  // new start
		loop_def.size()+1, // new end
		loop_def.cut() - loop_def.start() + 2, //new cut
		0.0 //skip_rate
	);
	set_loop( short_loop );

	using namespace fragment;
	FragSetOP short_frags = fragset->empty_clone(); // a fragset of same type should be able to handle everything
	FrameList loop_frames;
	fragset->region_simple( loop_def.start(), loop_def.stop(), loop_frames );
	for ( FrameList::const_iterator it = loop_frames.begin(),
			eit = loop_frames.end(); it!=eit; ++it ) {
		FrameOP short_frame = (*it)->clone_with_frags();
		short_frame->shift_to( short_frame->start() - loop_offset );
		short_frags->add( short_frame );
		tr.Trace << "ShortLoopClosure: short_frame->start" << short_frame->start() << std::endl;
		tr.Trace << "add Frame for ShortLoopClosure: "; short_frame->show_header( tr.Trace );
	}
	set_fragset( short_frags );

	init();
}

bool ShortLoopClosure::apply( Pose const& pose ) {
	mc_->set_update_boinc(false); // Dont send poses to observer cos they're truncated and the energies are meaningless.
	//copy loop segment into special purpose short_pose
	Pose short_pose;
	Size const short_size(  orig_loop_.size() + 2 );
	std::string sequence = pose.sequence();
	core::pose::make_pose_from_sequence(
		short_pose,
		sequence.substr( orig_loop_.start() - 2, orig_loop_.size() + 2 ),
		*( pose.residue_type_set_for_pose( chemical::CENTROID_t ))
	);
	short_pose.copy_segment( short_size, pose, 1, orig_loop_.start() - 1 );

	// make one jump fold-tree that connects the two fixed padding residues on either end
	kinematics::FoldTree f( short_size );
	f.new_jump( 1, short_size, loop().cut() );
	short_pose.fold_tree( f );
	tr.Debug << "FoldTree for short-pose segment in ShortLoopClosure "<< f << std::endl;
	// do the loop closure thing on this construct
	return Parent::apply( short_pose );
}

void ShortLoopClosure::catch_fragment( Pose const& short_pose ) {
	fragment::FrameOP aFrame = closure_fragments();
	aFrame->shift_to( 2 );
	aFrame->steal( short_pose );
	aFrame->shift_to( orig_loop_.start() );
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
