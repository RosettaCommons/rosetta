// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockRbTrialMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/MonteCarlo.hh>

// Package headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


// options
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

static basic::Tracer TR("protocols.moves.symmetry.SymFoldandDockRbTrialMover");

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover() :
	Mover( "SymFoldandDockRbTrialMover" ),
	smooth_move_(false), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn
) :
	Mover( "SymFoldandDockRbTrialMover" ),
	scorefxn_(scorefxn), smooth_move_(false), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move
) :
	Mover( "SymFoldandDockRbTrialMover" ),
	scorefxn_(scorefxn), smooth_move_(smooth_move), rot_mag_(8.0), trans_mag_(3.0), rigid_body_cycles_(50), mc_filter_(true)
{}

SymFoldandDockRbTrialMover::SymFoldandDockRbTrialMover(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool smooth_move,
	core::Real rot_mag,
	core::Real trans_mag
) :
	Mover( "SymFoldandDockRbTrialMover" ),
	scorefxn_( scorefxn), smooth_move_(smooth_move), rot_mag_(rot_mag), trans_mag_(trans_mag)
{}

void
SymFoldandDockRbTrialMover::apply( core::pose::Pose & pose )
{
	protocols::simple_moves::symmetry::SetupForSymmetryMover setup;
	setup.apply( pose );

	using namespace core::conformation::symmetry;
	using namespace basic::options;
	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	TR.Debug << "Rb move applied..." << std::endl;

	core::Real trans_mag_smooth = 0.1;
	core::Real rot_mag_smooth = 1.0;

	// Docking options
	if ( smooth_move_ ) {
		if ( option[ OptionKeys::fold_and_dock::trans_mag_smooth ].user() ) {
			trans_mag_smooth = option[ OptionKeys::fold_and_dock::trans_mag_smooth ];
		}
		if ( option[ OptionKeys::fold_and_dock::rot_mag_smooth ].user() ) {
			rot_mag_smooth = option[ OptionKeys::fold_and_dock::rot_mag_smooth ];
		}
	}

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rb_rot_magnitude ].user() )
	{
		rot_mag_ = ( option[ OptionKeys::fold_and_dock::rb_rot_magnitude ] );
	}

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rb_trans_magnitude ].user() )
	{
		trans_mag_ = ( option[ OptionKeys::fold_and_dock::rb_trans_magnitude ] );
	}

	core::Real rot_mag_trial = smooth_move_ ? rot_mag_smooth : rot_mag_;
	core::Real trans_mag_trial = smooth_move_ ? trans_mag_smooth : trans_mag_;

	// overrriding constructor or default values
	if ( option[ OptionKeys::fold_and_dock::rigid_body_cycles ].user() )
	{
		rigid_body_cycles_ = ( option[ OptionKeys::fold_and_dock::rigid_body_cycles ] );
	}

	if ( option[ OptionKeys::fold_and_dock::rigid_body_disable_mc ].user() )
	{
		mc_filter_ = false;
	}

	// Setup Monte Carlo object
	protocols::moves::MonteCarloOP monteCarlo_ = new protocols::moves::MonteCarlo(pose, *scorefxn_, 2.0 );

	// set up mover for docking
	protocols::rigid::RigidBodyDofSeqPerturbMover rb_perturb =
		protocols::rigid::RigidBodyDofSeqPerturbMover( dofs , rot_mag_trial, trans_mag_trial );

	if ( option[ OptionKeys::fold_and_dock::rotate_anchor_to_x ].user() ) {
		core::pose::symmetry::rotate_anchor_to_x_axis( pose );
	}

	for ( Size i = 1; i <= rigid_body_cycles_; ++i ) {
		rb_perturb.apply( pose );
		if ( mc_filter_ ) monteCarlo_->boltzmann( pose );
	}
	//monteCarlo_->recover_low(pose);
}

std::string
SymFoldandDockRbTrialMover::get_name() const {
	return "SymFoldandDockRbTrialMover";
}


} // symmetric_docking
} // protocols
