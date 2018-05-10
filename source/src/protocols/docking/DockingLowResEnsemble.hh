// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/DockinLowResEnsemble
/// @brief low resolution mode for ensemble docking
/// @details
/// @author Daisuke Kuroda


#ifndef INCLUDED_protocols_docking_DockingLowRes_Ensemble_hh
#define INCLUDED_protocols_docking_DockingLowRes_Ensemble_hh

#include <protocols/docking/types.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/DockingLowResEnsemble.fwd.hh>
#include <protocols/docking/DockingEnsemble.fwd.hh>

// Package headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/docking/ConformerSwitchMover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

#include <string>

// option key includes


#include <utility/vector1.hh>
#include <iostream>

namespace protocols {
namespace docking {

class DockingLowResEnsemble : public DockingLowRes
{
	//typedef core::Real Real;
public:

	/// @brief Default constructor
	DockingLowResEnsemble();

	// destructor
	~DockingLowResEnsemble() override;

	/// @brief Constructor with two arguments.  The first is scorefunction to be used for docking, the second is the
	///  DockJumps.
	DockingLowResEnsemble(
		core::scoring::ScoreFunctionCOP scorefxn,
		DockJumps const movable_jumps
	);

	void ensemble_defaults();

	protocols::moves::MoverOP clone() const override;
	// accessors
	core::Size get_current_conformer_ensemble1() override;
	core::Size get_current_conformer_ensemble2() override;

	void show( std::ostream & out=std::cout ) const override;

	// option setters
	void set_ensemble1( DockingEnsembleOP ensemble1 );
	void set_ensemble2( DockingEnsembleOP ensemble2 );

	void lowres_inner_cycle( core::pose::Pose & pose ) override;

	void adaptive_inner_cycle( core::pose::Pose & pose );
	void try_all_conformers( core::pose::Pose & pose );

	void initialize_movers();

	void calc_ensemble_accept_rate_();

	void reset_ensemble_trial_counters();


protected:
	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply

	void finalize_setup( core::pose::Pose & pose) override; // Comment out by DK

private:
	protocols::docking::ConformerSwitchMoverOP ensemble1_mover_;
	protocols::docking::ConformerSwitchMoverOP ensemble2_mover_;
	core::Size current_outer_cycle_;
	core::Real accept_rate_ens1_;
	core::Real num_trials_ens1_;
	core::Real num_accepted_ens1_;
	core::Real accept_rate_ens2_;
	core::Real num_trials_ens2_;
	core::Real num_accepted_ens2_;

	core::Real rb_trial_wt_;
	core::Real ens1_trial_wt_;
	core::Real ens2_trial_wt_;

	// Added by SSRB
	protocols::moves::SequenceMoverOP ensemble1_seq_;
	protocols::moves::SequenceMoverOP ensemble2_seq_;
	protocols::moves::TrialMoverOP ensemble1_trial_;
	protocols::moves::TrialMoverOP ensemble2_trial_;
	core::Size current_conformer_ensemble1_;
	core::Size current_conformer_ensemble2_;

	void ensemble1_trial( core::pose::Pose & pose );
	void ensemble2_trial( core::pose::Pose & pose );

};

} // docking
} // protocols

#endif
