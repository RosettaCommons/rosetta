// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/make_loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loop_build/LoopBuildMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/grafting/simple_movers/DeleteRegionMoverCreator.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <protocols/indel/IndelOptimizationMover.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingPrepackProtocol.hh>

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>

// C++ headers
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace indel {

static THREAD_LOCAL basic::Tracer TR("IndelOptimizationMover");

void
IndelOptimizationMover::apply(
	core::pose::Pose & pose
) {
	using namespace core;
	using namespace scoring;
	using namespace constraints;
	using namespace func;

	ScoreFunctionOP score_fxn = get_score_function();
	add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	add_fa_constraints_from_cmdline_to_pose(pose);

	pose.conformation().delete_residue_range_slow( start_res_, end_res_ );
	pose.conformation().detect_disulfides();
	pose.conformation().declare_chemical_bond( start_res_-1, "C", end_res_, "N" );
	TR << pose.fold_tree();

	// Add constraint
	/*if ( resnum > 1 ) {
	TR << "adding AtomPairConstraint between residues " << resnum-1 << " and " << resnum << std::endl;

	HarmonicFuncOP harm_func( new HarmonicFunc( 1.33, .02 ) );
	AtomID aidC( pose.residue( resnum-1 ).atom_index("C"), resnum-1 );
	AtomID aidN( pose.residue( resnum ).atom_index("N"), resnum );
	ConstraintCOP atompair( ConstraintOP( new AtomPairConstraint( aidC, aidN, harm_func ) ) );
	pose.add_constraint( atompair );
	score_fxn->set_weight( atom_pair_constraint, 1 );
	}*/

	// Loop model the remaining segment; range size before and
	Size loop_start = ( start_res_ < 2 + loop_length_ ) ? 2 : start_res_ - loop_length_;
	Size loop_end   = ( end_res_ + loop_length_ > pose.size() ) ? pose.size() : end_res_ + loop_length_;

	TR << "loop is from " << loop_start << " to " << loop_end << std::endl;

	bool there_is_a_jump = ( pose.conformation().num_chains() > 1 );

	std::string partners( &pose.pdb_info()->chain( loop_start ) );
	if ( there_is_a_jump ) {
		using namespace basic::options;
		if ( option[ OptionKeys::docking::partners ].user() ) {
			partners = option[ OptionKeys::docking::partners ].value();
		}
		/*.append( "_" );
		std::set< char > other_chains;
		//utility::vector1< Size > chain_endings = pose.conformation().chain_endings();
		for ( Size i = 1; i <= pose.size(); ++i ) {

		TR << "Residue " << i << " chain " << pose.pdb_info()->chain( i ) << std::endl;
		if ( pose.pdb_info()->chain( i ) != pose.pdb_info()->chain( loop_start ) ) {
		other_chains.insert( pose.pdb_info()->chain( i ) );
		TR << "Residue " << i << " chain " << pose.pdb_info()->chain( i ) << std::endl;
		}
		}
		for ( std::set< char >::iterator it = other_chains.begin(),
		end = other_chains.end(); it != end; ++it ) {
		partners.append( &(*it) );
		}
		*/
		TR << "Partners will be " << partners << std::endl;
	}

	Size cutpoint = Size( ( loop_start + loop_end ) / 2 );
	if ( cutpoint == loop_start ) ++cutpoint;
	if ( cutpoint == loop_end   ) --cutpoint;
	//protocols::loops::Loop loop( loop_start, loop_end, cutpoint );
	//TR << "i.e. " << loop << std::endl;
	loops::LoopsOP loops( new protocols::loops::Loops );
	loops->push_back( loop_start, loop_end, cutpoint );//loop );
	TR << (*loops) << std::endl;

	// Add constraints to native except loop
	pose::PoseOP native_pose( new core::pose::Pose( pose ) );
	native_pose->conformation().delete_residue_range_slow(loop_start, loop_end);
	native_pose->conformation().detect_disulfides();

	relax::AtomCoordinateCstMover coord_cst;
	coord_cst.set_refstruct( native_pose );
	coord_cst.apply( pose );

	// fragment initialization (not necessary in this case)
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( frag_files_ ) {
		loops::read_loop_fragments( frag_libs );
	}

	//setup of looprelax_mover
	comparative_modeling::LoopRelaxMover looprelax_mover;
	looprelax_mover.loops( loops );
	TR << ( *looprelax_mover.get_loops() ) << std::endl;
	looprelax_mover.frag_libs( frag_libs );
	looprelax_mover.relax( relax_ );
	looprelax_mover.refine( refine_ );
	looprelax_mover.remodel( remodel_ );
	looprelax_mover.intermedrelax( intermedrelax_ );

	loop_build::LoopBuildMoverOP loopbuild_mover( new loop_build::LoopBuildMover(looprelax_mover) );
	pose::Pose pose_saved( pose );
	loopbuild_mover->apply( pose );

	TR << "Finished first structure. Score is: " << ( *score_fxn )( pose ) << std::endl;

	jd2::JobOP job( jd2::JobDistributor::get_instance()->current_job() );

	if ( dump_initial_results_ ) {
		std::stringstream fn;
		fn << "pose_" << job->nstruct_index() << "_firstpass_1.pdb";
		pose.dump_pdb( fn.str() );
	}
	// If there is a jump, we may want to reoptimize the jump!
	// because our new structure is different

	// So grab the chain of the loop and build a docking protocol between that chain and the rest of the pose
	// Importantly, I think that we should actually keep 10 structures from the last operation
	// and dock each of them and return the best. Arbitrary number for now.
	utility::vector1< pose::Pose > poses;
	if ( there_is_a_jump ) {
		poses.push_back( pose );
		// Do the others
		for ( Size i = 2; i <= num_to_dock_; ++i ) {
			TR << "Okay, doing another structure ( " << i << " / " << num_to_dock_ << ")" << std::endl;
			pose::Pose temp_pose( pose_saved );
			loopbuild_mover->apply( temp_pose );
			poses.push_back( temp_pose );

			TR << "Finished another structure. Score is: " << ( *score_fxn )( temp_pose ) << std::endl;

			if ( dump_initial_results_ ) {
				std::stringstream fn;
				fn << "pose_" << job->nstruct_index() << "_firstpass_" << i << ".pdb";
				temp_pose.dump_pdb( fn.str() );
			}
		}

		for ( Size i = 1; i <= num_to_dock_; ++i ) {
			TR << "Score for pose " << i << " is " << ( *score_fxn )( poses[ i ] ) << std::endl;
		}

		// What jumps should be moveable? The ones that lead to the Special Chain
		// Therefore, the jumps for which the upstream residue or downstream residue are Special Chain
		utility::vector1< Size > moveable_jumps;
		Size n_jumps = pose.fold_tree().num_jump();
		for ( Size i = 1; i <= n_jumps; ++i ) {
			Size resi = pose.fold_tree().upstream_jump_residue( i );
			if ( pose.pdb_info()->chain( resi ) == pose.pdb_info()->chain( loop_start ) ) {
				moveable_jumps.push_back( i );
				continue;
			}
			resi = pose.fold_tree().downstream_jump_residue( i );
			if ( pose.pdb_info()->chain( resi ) == pose.pdb_info()->chain( loop_start ) ) {
				moveable_jumps.push_back( i );
				continue;
			}
		}

		docking::DockingPrepackProtocolOP ppk( new docking::DockingPrepackProtocol() );
		ppk->set_partners( partners );
		docking::DockingProtocolOP dock( new docking::DockingProtocol( moveable_jumps, false, true ) );
		dock->set_partners( partners );

		// Create a docking prepack protocol and a docking protocol

		// Now apply all our docking stuff
		Real min_score = 100000;
		Size best_i = 1;
		for ( Size i = 1; i <= poses.size(); ++i ) {
			TR << "Prepacking and docking result " << i << " / " << poses.size() << std::endl;
			ppk->apply( poses[ i ] );
			dock->apply( poses[ i ] );
			Real score = ( *score_fxn )( poses[ i ] );
			if ( score < min_score ) {
				best_i = i;
			}
		}

		pose = poses[ best_i ];
	}

	return;
}

}
}
