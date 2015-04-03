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
/// @author

#include <devel/dna/ProteinDNA_Relax.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/Tracer.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/rigid/rigid_body_moves.hh>
#include <protocols/moves/RepeatMover.hh>

// ObjexxFCL Headers

// Utility Headers

// c++ headers
#include <iostream>
#include <string>
#include <vector>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using utility::vector1;

namespace devel {
namespace dna {

using namespace core;
using namespace conformation;
using namespace chemical;
// file-scope -- prob bad
static thread_local basic::Tracer tt( "devel.dna.ProteinDNA_Relax", basic::t_info );


////////////////////////////////////////////////////////////////////////////////////////////////////

protocols::moves::TrialMoverOP
setup_MCM_trial(
	protocols::moves::MoverOP perturb,
	protocols::moves::MoverOP pack,
	protocols::moves::MoverOP minimize,
	protocols::moves::MonteCarloOP mc
)
{
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( perturb  );
	seq->add_mover( pack     );
	seq->add_mover( minimize );

	return protocols::moves::TrialMoverOP( new TrialMover( seq, mc ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void
RB_Mover::apply( core::pose::Pose & pose )
{
	moved_jump_ = protocols::rigid::gaussian_jump_move( pose, *mm_, trans_mag_, rot_mag_ );
}

std::string
RB_Mover::get_name() const {
	return "RB_Mover";
}


////////////////////////////////////////////////////////////////////////////////////////////////////

ProteinDNA_Relax::ProteinDNA_Relax( ScoreFunction const & scorefxn_in, Size const moving_jump_in ):
	protocols::moves::Mover( "ProteinDNA_Relax" ),
	scorefxn_( scorefxn_in.clone() ),
	moving_jumps_( 1, moving_jump_in )
{
	apply_default_settings();
}

ProteinDNA_Relax::ProteinDNA_Relax( ScoreFunction const & scorefxn_in, utility::vector1< Size > const & moving_jumps ):
	protocols::moves::Mover( "ProteinDNA_Relax" ),
	scorefxn_( scorefxn_in.clone() ),
	moving_jumps_( moving_jumps )
{
	apply_default_settings();
}

void
ProteinDNA_Relax::apply( pose::Pose & pose )
{
	using namespace protocols::moves;

	Size const nres( pose.total_residue() );

	// setup task and move maps
	kinematics::MoveMapOP mm( new kinematics::MoveMap ); // default state is DOF FIXED
	pack::task::PackerTaskOP pack_task( pack::task::TaskFactory::create_packer_task( pose ));

	pack_task->initialize_from_command_line();
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			pack_task->nonconst_residue_task( i ).restrict_to_repacking();
			mm->set_chi( i, true );

		} else {
			pack_task->nonconst_residue_task( i ).restrict_to_repacking();
			mm->set_chi( i, true );

		}
	}
	pack_task->or_include_current( true );

	for ( Size i=1; i<= moving_jumps_.size(); ++i ) {
		mm->set_jump( moving_jumps_[i], true );
	}

	// use more rotamers for rotamer trials, since its faster
	pack::task::PackerTaskOP rottrial_task( pack_task->clone() );

	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			rottrial_task->nonconst_residue_task(i).or_ex1( true );
			rottrial_task->nonconst_residue_task(i).or_ex2( true );
			rottrial_task->nonconst_residue_task(i).or_ex1aro( true );
			rottrial_task->nonconst_residue_task(i).or_ex1aro_sample_level( pack::task::EX_THREE_THIRD_STEP_STDDEVS );
			rottrial_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		}
	}

	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn_, 0.8 ) );

	// this is really annoying -- ask Andrew how to fix:
	(*scorefxn_)(pose);

	//////////////////
	// setup movers:

	// packmover
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( scorefxn_, pack_task, 25 ) );

	// min mover
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn_, "dfpmin", min_tol_, true ) );

	// rb mover
	RB_MoverOP rb_mover( new RB_Mover( mm, trans_mag_, rot_mag_ ) );

	// rotamer trials w/ energycut
	protocols::simple_moves::EnergyCutRotamerTrialsMoverOP rottrial_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn_, *rottrial_task, mc, energycut_ ) );

	// trials:
	TrialMoverOP rb_min_trial = setup_MCM_trial( rb_mover, rottrial_mover, min_mover, mc );

	TrialMoverOP pack_trial( new TrialMover( pack_mover, mc ) );

	TrialMoverOP min_trial( new TrialMover( min_mover, mc ) );

	// now the inner cycles mover
	RepeatMoverOP rb_min_cycle( new RepeatMover( rb_min_trial, inner_cycles_ ) );

	// ramp up repulsive, dunbrack(?)
	//
	// this guy to the (ramp_cycles-1) power should be 1.0 / initial_weight
	core::Real const ramping_multiplier( ( ramping_cycles_ < 1 ) ?
																			 1.0 : std::exp( std::log( 1.0 / ramping_initial_weight_ ) / ( ramping_cycles_)));
	core::Real ramping_weight_factor( ramping_initial_weight_ / ramping_multiplier );
	core::Real const final_fa_rep_weight( (*scorefxn_)[ scoring::fa_rep ] );

	for ( Size m1=1; m1<= outer_cycles_; ++m1 ) {

		if ( ramping_cycles_ && m1 <= ramping_cycles_+1 ) {
			// ramp the weights
			ramping_weight_factor *= ramping_multiplier;
			scorefxn_->set_weight( scoring::fa_rep, final_fa_rep_weight * ramping_weight_factor );
			mc->score_function( *scorefxn_ );

		} else {
			assert( std::abs( ramping_weight_factor - 1.0 ) < 1e-3 );
		}

		pack_trial->apply( pose );

		min_trial->apply( pose );

		rb_min_cycle->apply( pose );

		mc->recover_low( pose );

		mc->show_counters();
	}

	pack_trial->apply( pose );

	min_trial->apply( pose );

	mc->recover_low( pose );


}

std::string
ProteinDNA_Relax::get_name() const {
	return "ProteinDNA_Relax";
}


} // namespace dna
} // namespace devel
