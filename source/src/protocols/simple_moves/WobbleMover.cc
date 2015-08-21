// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Oliver Lange

// Unit Headers
#include <protocols/simple_moves/WobbleMover.hh>

// Package Headers
#include <protocols/simple_moves/GunnCost.hh>

// Project Headers
#include <core/fragment/FragSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>

// Utility headers
#include <numeric/random/random.hh>

#include <utility/vector1.hh>


// C++ headers

static thread_local basic::Tracer TR( "protocol.abinitio.WobbleMover" );

namespace protocols {
namespace simple_moves {


using namespace core;
using namespace fragment;

WobbleMover::WobbleMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap,
	protocols::simple_moves::FragmentCostOP cost )
:
	protocols::simple_moves::ClassicFragmentMover( fragset, movemap, "WobbleMover" ), // explicit initialization of virtual base class required
	protocols::simple_moves::SmoothFragmentMover( fragset, movemap, cost, "WobbleMover")
{
	set_defaults();
}

WobbleMover::WobbleMover(
	core::fragment::FragSetCOP fragset,
	core::kinematics::MoveMapCOP movemap
) :
	protocols::simple_moves::ClassicFragmentMover( fragset, movemap, "WobbleMover" ), // explicit initialization of virtual base class required
	protocols::simple_moves::SmoothFragmentMover( fragset, movemap, FragmentCostOP( new protocols::simple_moves::GunnCost ), "WobbleMover")
{
	set_defaults();
}

/// Copy constructor disabled until the virtual base class
/// ClassicFragmentMover's construction is sorted out
/*WobbleMover::WobbleMover( WobbleMover const & src ) :
Parent( src ),
buffer_length_( src.buffer_length_ ),
forward_threshold_( src.forward_threshold_ ),
backward_threshold_( src.backward_threshold_ )
{}*/

WobbleMover::~WobbleMover()
{}

std::string
WobbleMover::get_name() const {
	return "WobbleMover";
}

bool WobbleMover::ccd_closure(
	pose::Pose & pose,
	protocols::loops::Loops const & loops,
	kinematics::MoveMap const & mm ) const
{
	// There is only one loop.
	protocols::loops::Loops::const_iterator it = loops.begin();
	protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover(
		*it, kinematics::MoveMapCOP( kinematics::MoveMapOP( new kinematics::MoveMap( mm ) ) ) );
	ccd_loop_closure_mover.apply( pose );

	return ( ccd_loop_closure_mover.deviation() < forward_threshold_ );
}

/// @brief make a wobble move ( smooth move + ccd loop closure )
/// @detail NEEDS MORE WORK IF MULTIPLE CHAINS OR JUMPS ARE USED: total_residue!=end of chain
bool WobbleMover::apply_fragment (
	Frame const& frame,
	Size frag_num,
	kinematics::MoveMap const& movemap,
	pose::Pose &pose ) const
{

	// define a loop with cut-point, such that fragment is on one side until cut-point and otherside has buffer_ residues.
	protocols::loops::Loops frag_loop;


	// if fragment is exactly at start or end of pose it is just inserted without ccd
	bool use_ccd ( ! ( frame.start() == 1 || frame.end() >= pose.total_residue()-1 ) );
	kinematics::FoldTree original_tree;
	TR.Debug << " start : " << frame.start() <<" buffer " << frame.start()-buffer_length_ << " end: " << frame.end() << " total res: " << pose.total_residue() << std::endl;
	if ( use_ccd ) {
		// prepare foldtree for ccd

		bool cut_Cterm; // controls whether cutpoint is at C-terminal side of fragment
		if ( frame.end() + buffer_length_ >= pose.total_residue() ) { //close to end of pose: cut at Nterm
			cut_Cterm = false;
		} else if ( frame.start() <= buffer_length_ ) { //close to start of pose: cut at Cterm
			cut_Cterm = true;
		} else { //otherwise random direction
			cut_Cterm = numeric::random::rg().uniform() >= 0.5 ;
		};

		if ( cut_Cterm ) {
			frag_loop.add_loop( frame.start()-1, frame.end()+buffer_length_, frame.end() );
		} else {
			frag_loop.add_loop( frame.start()-buffer_length_, frame.end()+1, frame.start()-1 );
		};

		TR.Debug << "loop definitio for wobble move" << std::endl << frag_loop << std::endl;

		//set fold-tree ( keep original for later )
		original_tree = pose.fold_tree();
		TR.Debug << "original " << original_tree << std::endl;
		kinematics::FoldTree fold_tree;
		protocols::loops::fold_tree_from_loops( pose, frag_loop, fold_tree, true );
		TR.Debug << "for loops " << fold_tree << std::endl;
		pose.fold_tree( fold_tree );
		loops::add_cutpoint_variants( pose );

	} // ccd preparation

	//apply fragment
	bool success = frame.apply( movemap, frag_num, pose );
	// pose.dump_pdb( "pre_ccd.pdb" );

	//close loop with ccd
	if ( success && use_ccd ) {
		pose::Pose orig_pose = pose;
		success = ccd_closure( pose, frag_loop, movemap );
		if ( !success ) pose = orig_pose;
	}

	//clean up
	if ( use_ccd ) {
		loops::remove_cutpoint_variants( pose );
		pose.fold_tree( original_tree );
	}

	return success;
}

} //abinitio
} //protocols
