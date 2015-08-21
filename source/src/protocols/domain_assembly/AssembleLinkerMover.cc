// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/AssembleLinkerMover.cc
/// @brief  A mover for assembling domains
/// @author James Thompson

#include <protocols/domain_assembly/AssembleLinkerMover.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>


#include <core/chemical/ChemicalManager.fwd.hh>


#include <protocols/moves/Mover.hh>

//#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes


#include <core/util/SwitchResidueTypeSet.hh>


namespace protocols {
namespace domain_assembly {

void AssembleLinkerMover::apply( core::pose::Pose & pose ) {
	using core::Real;
	using std::string;

	char const first_chain( pose.pdb_info()->chain(1) );
	// rename second chain
	Size breakpoint(0);
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		char const chain_ii( pose.pdb_info()->chain(ii) );
		if ( !breakpoint && chain_ii != first_chain ) {
			breakpoint = ii;
		}
		pose.pdb_info()->number(ii,ii);
	}
	// make a single-chain pdb
	pose.pdb_info()->set_chains(first_chain);
	pose.conformation().reset_chain_endings();

	// new fold tree
	using core::kinematics::FoldTree;
	FoldTree new_fold_tree(pose.total_residue());
	pose.fold_tree(new_fold_tree);

	// remodel loops
	using namespace protocols::loops;
	using namespace protocols::comparative_modeling;
	Size const loop_start( breakpoint - min_loop_size_ );
	Size const loop_stop ( breakpoint + min_loop_size_ );

	core::util::switch_to_residue_type_set(
		pose, core::chemical::CENTROID
	);

	Loop loop( loop_start, loop_stop, 0, 0, false );
	LoopsOP loops( new Loops() );
	loops->add_loop(loop);
	loop_mover::LoopMoverOP loop_mover = LoopMoverFactory::get_instance()->create_loop_mover(
		loop_mover_name_, loops
	);
	for ( Size ii = 1; ii <= frag_libs_.size(); ++ii ) {
		loop_mover->add_fragments( frag_libs_[ii] );
	}
	loop_mover->apply(pose);
	Size const max_tries( 10 ); // make this an option?
	bool loops_closed( false );
	for ( Size ii = 1; (ii <= max_tries) && !loops_closed; ++ii ) {
		LoopsOP loops = pick_loops_chainbreak( pose, min_loop_size_ );
		loops_closed = ( loops->size() == 0 );
		if ( loops_closed ) {
			loop_mover = LoopMoverFactory::get_instance()->create_loop_mover( loop_mover_name_, loops );
			loop_mover->apply( pose );
		}
	}

	core::util::switch_to_residue_type_set(
		pose, core::chemical::FA_STANDARD
	);
} // apply


}//namespace domain_assembly
}//namespace protocols

