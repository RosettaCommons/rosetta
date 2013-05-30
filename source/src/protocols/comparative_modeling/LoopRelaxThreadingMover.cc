// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Mike Tyka

// include these first for building on Visual Studio

#include <protocols/loops/loops_main.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/comparative_modeling/LoopRelaxThreadingMover.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
// AUTO-REMOVED #include <protocols/electron_density/util.hh>

#include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <utility/vector1.hh>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
// AUTO-REMOVED #include <protocols/jd2/MultiThreadingJob.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>
// AUTO-REMOVED #include <protocols/comparative_modeling/MultiThreadingMover.hh>
#include <protocols/simple_moves/RepulsiveOnlyMover.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace comparative_modeling {

void LoopRelaxThreadingMover::setup() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// initialize fragments
	loops::read_loop_fragments( frag_libs_ );

	// Read parameters
	max_loop_rebuild     = option[ OptionKeys::cm::max_loop_rebuild ]();
	loop_rebuild_filter  = option[ OptionKeys::cm::loop_rebuild_filter ]();
	remodel              = option[ OptionKeys::loops::remodel ]();
	relax                = option[ OptionKeys::loops::relax ]();
}

void LoopRelaxThreadingMover::apply( core::pose::Pose & pose ) {
	using namespace protocols::loops;
	using namespace protocols::jd2;

	// apply a mover which calculates only repulsive energy on designate residues
	protocols::simple_moves::RepulsiveOnlyMover replonly;
	replonly.apply( pose );

	using core::Size;
	basic::Tracer tr("protocols.threading");
	// looprelax
	protocols::comparative_modeling::ThreadingJobCOP job = dynamic_cast< protocols::comparative_modeling::ThreadingJob const * >(
		JobDistributor::get_instance()->current_job()->inner_job().get()
	);
	if ( !job ) {
		utility_exit_with_message("ERROR: You must use the ThreadingJobInputter with the LoopRelaxThreadingMover - did you forget the -in:file:template_pdb option?");
	}

	LoopsOP my_loops = new Loops( job->loops( pose.total_residue() ) );
	my_loops->choose_cutpoints( pose );
	tr.Info << "loops to be rebuilt are: " << std::endl;
	tr.Info << my_loops << std::endl;

	// Add any ligands specified
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// setup for symmetry
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
		protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
		pre_mover.apply( pose );
	}

	// Skip modelling if we're just producing starting models
	if ( option[ OptionKeys::cm::start_models_only ]() ) return;

	// setup for density
	if ( option[ edensity::mapfile ].user() ) {
		protocols::electron_density::SetupForDensityScoringMover pre_mover;
		pre_mover.mask( *my_loops );
		pre_mover.apply( pose );
	}

	LoopRelaxMoverOP lr_mover( new LoopRelaxMover );
	lr_mover->frag_libs( frag_libs_ );
	lr_mover->loops( my_loops );
	lr_mover->relax( relax );
	lr_mover->remodel( remodel );
	lr_mover->cmd_line_csts( true );
	lr_mover->rebuild_filter( loop_rebuild_filter );
	lr_mover->n_rebuild_tries( max_loop_rebuild );
	lr_mover->copy_sidechains( true );
	lr_mover->set_current_tag( get_current_tag() );
	lr_mover->apply( pose );

	// more loop rebuilding (if necessary!)
	Size const min_loop_size( option [ OptionKeys::cm::min_loop_size ]() );
	if ( option[ OptionKeys::cm::loop_close_level ]() == 2 ) {
		// if loops aren't closed here, try to figure out the loops again
		using namespace protocols::comparative_modeling;
		LoopsOP temp_loops = pick_loops_chainbreak( pose, min_loop_size );
		lr_mover->loops( temp_loops );
		lr_mover->apply( pose );
	} else if ( option[ OptionKeys::cm::loop_close_level ]() == 3 ) {
		using namespace protocols::comparative_modeling;
		Size const max_tries( 10 ); // make this an option?
		bool loops_closed( false );
		for ( Size ii = 1; (ii <= max_tries) && !loops_closed; ++ii ) {
			LoopsOP temp_loops = pick_loops_chainbreak( pose, min_loop_size );
			loops_closed = ( temp_loops->size() == 0 );
			if ( !loops_closed ) {
				lr_mover->loops(temp_loops);
				lr_mover->apply(pose);
				lr_mover->relax("no");
			}
		}

		lr_mover->remodel("no");
		lr_mover->relax(relax);
		lr_mover->apply(pose);
	}

	// recover side chains if requested by user ...
	if ( option[ cm::recover_side_chains ]() ) {
		using namespace core::pose;
		using namespace core::pack;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::sequence;
		using namespace protocols::comparative_modeling;

		Pose template_pose( *job->get_pose() );

		SequenceAlignment aln = alignment_from_pose(pose);
		SequenceMapping map   = aln.sequence_mapping(1,2);
		utility::vector1< bool > residues_to_repack( pose.total_residue(), true );

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			//if ( map[ii] == 0 ) residues_to_repack[ii] = true;
			residues_to_repack[ii] = true;
		}
		StealSideChainsMover sc_mover( template_pose, map );
		sc_mover.apply( pose );
		task::PackerTaskOP task
			= task::TaskFactory::create_packer_task( pose );
		task->initialize_from_command_line();
		task->restrict_to_repacking();
		task->restrict_to_residues(residues_to_repack);
		ScoreFunctionOP scorefxn(
			getScoreFunction()
		);
		pack_rotamers( pose, *scorefxn, task );
	} // recover_side_chains
} // apply

std::string
LoopRelaxThreadingMover::get_name() const {
	return "LoopRelaxThreadingMover";
}

} // namespace loops
} // namespace protocols
