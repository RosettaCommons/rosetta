// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/LoopJumpFoldCst.hh>

// Package Headers
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>


#include <basic/options/option.hh>


#include <protocols/loops/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/simple_moves/FragmentMover.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/util.hh>
#include <basic/Tracer.hh>


//// C++ headers
#include <fstream>


// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.general_abinitio", basic::t_info );

namespace protocols {
namespace abinitio {

using namespace core;

LoopJumpFoldCst::~LoopJumpFoldCst() = default;

std::string
LoopJumpFoldCst::get_name() const {
	return "LoopJumpFoldCst";
}

///////////////////////////////////////////////////////////////////////
/// @brief Select loop set at random using skip rate
void LoopJumpFoldCst::select_loops(
	loops::Loops& loops_out
) const
{
	loops_out.clear();
	int ntries = 0;
	while ( loops_out.size() == 0 && ntries++ < 50 ) {
		for ( auto const & loop : loops_ ) {
			if ( numeric::random::rg().uniform() >= loop.skip_rate() )  {
				loops_out.push_back( loop );
			}
		}
	}
	if ( loops_out.size() == 0 ) {
		loops_out = loops_;
	}
} // void LoopRebuild::select_loops

// helper to new_kinematics
bool LoopJumpFoldCst::parse_jump_def( KinematicControlOP current_kinematics, kinematics::MoveMapOP movemap ) {
	// get flexible jumps ( beta-sheet stuff etc. )
	jumping::JumpSample current_jumps;
	if ( jump_def_ ) {
		movemap->set_jump( true ); //careful that these don't get minimized!
		// initialize jumping
		Size attempts( 10 );
		do {
			current_jumps = jump_def_->create_jump_sample();
		} while ( !current_jumps.is_valid() && attempts-- );

		if ( !current_jumps.is_valid() ) {
			utility_exit_with_message( "not able to build valid fold-tree in JumpingFoldConstraints::setup_foldtree" );
		}

		tr.Debug << "current_jumps " << current_jumps << std::endl;

		// get fragments to sample these jumps
		core::fragment::FragSetOP jump_frags;

		if ( coordinate_constraint_weight_ > 0 ) { //don't allow jumps move in xyz-fixed area --- i.e., the stuff that should stay rather fixed
			// movemap = new kinematics::MoveMap( *movemap ); seems redundant
			// flexible_part.switch_movemap( *movemap, id::BB, true ); this line seems redundant here
		}
		jump_frags = jump_def_->generate_jump_frags( current_jumps, *movemap );
		simple_moves::ClassicFragmentMoverOP jump_mover( new protocols::simple_moves::ClassicFragmentMover( jump_frags, movemap ) );
		jump_mover->type( "JumpMoves" );
		jump_mover->set_check_ss( false ); // this doesn't make sense with jump fragments
		jump_mover->enable_end_bias_check( false ); //no sense for discontinuous fragments


		current_kinematics->set_sampling_fold_tree( current_jumps.fold_tree() );
		current_kinematics->set_jump_mover( jump_mover );
	} else {
		//take existing fold-tree
		current_kinematics->set_sampling_fold_tree( current_kinematics->final_fold_tree() );
	}

	return true;
}


//@brief create a new fold-tree and movemap --- a KinematicControl object
// a basic generalized protocol:
// loops are determined: if loops present only loops move, otherwise everything moves
// get jumps from JumpDef ( as in JumpingFoldCst, for instance beta-sheet jumps )
// combine jumps with additional jumps that keep the rigid portions together ( looprlx-type )
// set movemap to allow only sampling of flexible parts. call sampling protocol
KinematicControlOP LoopJumpFoldCst::new_kinematics( pose::Pose &pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// select loops
	loops::Loops loops;
	if ( loops_.size() > 0 ) {
		select_loops( loops );
	}
	tr.Debug << loops << std::endl;

	KinematicControlOP current_kinematics( nullptr );
	if ( loops.size() && coordinate_constraint_weight_ > 0.0 ) {
		current_kinematics = KinematicControlOP( new CoordinateConstraintKC( false /*ramp*/, coordinate_constraint_weight_ ) );
	} else {
		current_kinematics = KinematicControlOP( new KinematicControl );
	}

	loops::Loops mmloops;
	if ( option[ OptionKeys::loops::mm_loop_file ].user() ) {
		mmloops = loops::Loops( option[ OptionKeys::loops::mm_loop_file ]() );
	} else {
		mmloops = loops;
	}

	//figure out movemap!
	kinematics::MoveMapOP movemap( new kinematics::MoveMap );
	movemap->set_jump( true );
	if ( jump_def_ ) movemap->set_jump( true ); //careful that these don't get minimized!
	if ( mmloops.size() && coordinate_constraint_weight_ == 0.0 ) {
		mmloops.switch_movemap( *movemap, id::BB, true );
		loops.switch_movemap( *movemap, id::BB, true );
	} else {
		movemap->set_bb( true );
	}
	loops::Loops rigid_core( loops.invert( pose.size() ) );

	bool success( true );
	current_kinematics->set_final_fold_tree( pose.fold_tree() );
	success = parse_jump_def( current_kinematics, movemap );
	current_kinematics->set_movemap( movemap );

	// change fold-tree such that only loop parts can move
	if ( loops.size() ) {
		// subsequently add sufficiently many jumps such that all numbers are the same.
		success &= add_rigidity_jumps( rigid_core, current_kinematics );
		if ( !success ) {
			tr.Warning << "[WARNING] was not able to fix rigid regions with jumps...retry" << std::endl;
			return nullptr;
		}
	}
	pose.fold_tree( current_kinematics->sampling_fold_tree() );
	set_extended_torsions_and_idealize_loops( pose, loops );
	if ( rigid_core.size() ) {
		/*success = */add_coord_cst( rigid_core, pose );
	}
	return current_kinematics;
}

bool
LoopJumpFoldCst::add_coord_cst( loops::Loops const& rigid_core, core::pose::Pose& pose ) {
	if ( coordinate_constraint_weight_ != 0.0 ) {
		bool bWeightsIn = ( coordinate_constraint_weights_.size() >= pose.size() );
		loops::fix_with_coord_cst( rigid_core, pose, bCstAllAtom_, coordinate_constraint_weights_ );
		if ( dump_weights_file_ != "NO_DUMP" ) {
			utility::io::write_vector( dump_weights_file_, coordinate_constraint_weights_ );
			dump_weights_file_ = "NO_DUMP"; //do it only once
		}

		//add also constraints to the resolution switch...
		loops::fix_with_coord_cst( rigid_core, res_switch().init_pose(), true, coordinate_constraint_weights_ );
		if ( !bWeightsIn ) coordinate_constraint_weights_.clear(); //delete them again if they haven't been provided from the outside
	}
	return true;
}

bool
LoopJumpFoldCst::add_rigidity_jumps( loops::Loops const& rigid, KinematicControlOP current_kinematics ) {
	ObjexxFCL::FArray1D_float const& cut_probability( ss_def_->loop_fraction() );
	// add jumps to the fold-tree such that the rigid parts are all connected by jumps.
	// first find out which rigid regions ( which are nothing else than inverted loops )
	// are already connected by a jump.
	// we make a "visited" field. It will contain number 0 for regions that have no connection to rigid regions
	// and for each rigid-super-region it has a distinct number, e.g., 1 2 1 0 1 2 0, means region 1 3 with 5 and 2 with 6 are already connected
	kinematics::FoldTree const f( current_kinematics->sampling_fold_tree() );

	// list of jumps that connect rigid regions ( these jump_dofs should not be moved )
	jumping::JumpList rigid_jumps;
	// list of jumps that connect flexible regions ( these jump_dofs should be moved )
	jumping::JumpList flex_jumps;

	// go thru existing jumps and fill them into the rigid_jump/flex_jump lists. Moreover, figure out which rigid regions are already connected by jumps
	utility::vector1< Size > visited( rigid.size(), 0 );
	for ( Size jump_nr = 1; jump_nr <= f.num_jump(); jump_nr++ ) {
		Size const up( f.upstream_jump_residue( jump_nr ) );
		Size const down( f.downstream_jump_residue( jump_nr ));
		Size up_loop;
		Size down_loop;
		if ( (up_loop = rigid.loop_index_of_residue( up )) && (down_loop = rigid.loop_index_of_residue( down ) ) ) {
			if ( up_loop != down_loop
					&& ( !visited[ up_loop ] || !visited[ down_loop ] || ( visited[ up_loop ] != visited[ down_loop ] ) ) ) {
				rigid_jumps.push_back( jumping::Interval( up, down ) ); //good jump  --- keep it
				// decide upon 3 cases: both nodes unvisited, 1 node visited, both nodes visited
				// case0 : both new tag with new jump_nr
				Size visit_nr = jump_nr;
				// case1 : both visited--> replace all higher visit_nr by lower visit_nr
				if ( visited[ up_loop ] && visited[ down_loop ] ) {
					Size old_visit_nr = visited[ down_loop ]; //arbitrary choice
					visit_nr = visited[ up_loop ];
					for ( Size i=1; i<=visited.size(); i++ ) {
						if ( visited[ i ] == old_visit_nr ) visited[ i ] = visit_nr;
					}
				} else if ( visited[ up_loop ] || visited[ down_loop ] ) {
					// case2: one already visited the other is still zero and thus neutral to addition
					visit_nr = visited[ up_loop ] + visited[ down_loop ];
				}
				visited[ up_loop ] = visit_nr;
				visited[ down_loop ] = visit_nr;
			} // jump between different regions
		} else {
			flex_jumps.push_back( jumping::Interval( up, down ) ); // jump does not link rigid regions
		}
	} // for jump_nr

	//now we have a connection pattern based on the jumps already present.
	//a visited region and make it the root-reg
	Size root_reg = 0;
	for ( Size region = 1; region <= visited.size(); region++ ) {
		if ( visited[ region ] ) {
			root_reg = region;
			break;
		}
	}
	// if no rigid regions are yet connected, define one arbitrarily as the root-reg
	if ( root_reg == 0 ) {
		root_reg = 1;
		visited[ root_reg ] = 1;
	}

	for ( auto & rigid_jump : rigid_jumps ) {
		tr.Debug << "Fix_jumps: " << rigid_jump.start_<< " " << rigid_jump.end_ << std::endl;
	}
	tr.Debug << "now add more fix-jumps " << std::endl;

	loops::Loops::LoopList rigid_loops = rigid.loops(); // loops in sequence that correspond to the regions
	Size const anchor( static_cast< Size >( 0.5*(rigid_loops[ root_reg ].stop() - rigid_loops[ root_reg ].start()) ) + rigid_loops[ root_reg ].start() );
	for ( Size region = 1; region <= visited.size(); region++ ) {
		Size old_visited = visited[ region ];
		if ( visited[ region ] != visited[ root_reg ] ) {
			rigid_jumps.push_back ( jumping::Interval( anchor,
				rigid_loops[ region ].start() + static_cast< Size >( 0.5*( rigid_loops[ region ].stop()-rigid_loops[ region ].start() ) ) ) );
			visited[ region ] = visited[ root_reg ];
			if ( old_visited ) { //if we connected a cluster make sure to update all its nodes
				for ( Size i=1; i<=visited.size(); i++ ) {
					if ( visited[ i ] == old_visited ) visited[ i ] = visited[ root_reg ];
				}
			}
		}
	} // for region


	ObjexxFCL::FArray1D_float new_cut_prob( cut_probability );
	for ( auto const & it : rigid ) {
		for ( Size pos = it.start(); pos <= it.stop(); pos++ ) {
			new_cut_prob( pos ) = 0;
		}
	}

	Size total_njump( rigid_jumps.size() + flex_jumps.size() );
	ObjexxFCL::FArray2D_int jumps( 2, total_njump );
	Size ct = 1;
	for ( auto it = rigid_jumps.begin(), eit = rigid_jumps.end(); it !=eit; ++it, ++ct ) {
		jumps( 1, ct ) = it->start_;
		jumps( 2, ct ) = it->end_;
		tr.Debug << "Fix_jumps: " << it->start_<< " " << it->end_ << std::endl;
	}
	for ( auto it = flex_jumps.begin(), eit = flex_jumps.end(); it !=eit; ++it, ++ct ) {
		jumps( 1, ct ) = it->start_;
		jumps( 2, ct ) = it->end_;
		tr.Debug << "Flex_jumps: " << it->start_<< " " << it->end_ << std::endl;
	}
	for ( Size i = 1; i<= current_kinematics->sampling_fold_tree().nres(); i++ ) {
		tr.Trace << "cut_prob: " << i << " " << new_cut_prob( i ) << std::endl;
	}
	kinematics::FoldTree new_fold_tree;
	bool success = new_fold_tree.random_tree_from_jump_points( current_kinematics->sampling_fold_tree().nres(), total_njump, jumps, new_cut_prob );

	if ( success ) {
		if ( rigid_jumps.size() ) {
			//this reordering is essential if we use coord cst. since they are relative to the origin (root of fold-tree)
			// we want the root in the most "fixed" position, and then fix the coord-cst atoms to this root
			success = new_fold_tree.reorder( rigid_jumps.begin()->start_ ); //or try multiple residues?
			if ( !success ) {
				tr.Warning << "[WARNING] failed to reorder fold-tree() ... retry"<<std::endl;
				return success;
			}
		}
		new_fold_tree.put_jump_stubs_intra_residue();
		current_kinematics->set_sampling_fold_tree( new_fold_tree );
		kinematics::MoveMapOP  mm( new kinematics::MoveMap( current_kinematics->movemap() ) );
		mm->set_jump( false );
		//find flexible jump-nr in fold tree and update movemap accordingly
		for ( auto & flex_jump : flex_jumps ) {
			mm->set_jump( flex_jump.start_, flex_jump.end_, true );
		}
		current_kinematics->set_movemap( mm );
	}

	return success;
}

}
}
