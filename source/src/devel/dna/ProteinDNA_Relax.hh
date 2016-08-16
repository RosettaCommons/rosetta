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
/// @author

#ifndef INCLUDED_devel_dna_ProteinDNA_Relax_hh
#define INCLUDED_devel_dna_ProteinDNA_Relax_hh


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/TrialMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <string>

#include <utility/vector1.hh>


namespace devel {
namespace dna {


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

class RB_Mover : public protocols::moves::Mover {
public:
	typedef core::Real Real;

public:
	RB_Mover(
		core::kinematics::MoveMapOP const mm_in,
		Real const trans_mag_in,
		Real const rot_mag_in
	):
		protocols::moves::Mover( "RB_Mover" ),
		mm_( mm_in ),
		trans_mag_( trans_mag_in ),
		rot_mag_( rot_mag_in ),
		moved_jump_( 0 )
	{}


	void
	apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	Size
	moved_jump() const
	{
		return moved_jump_;
	}

private:
	core::kinematics::MoveMapOP const mm_;
	Real const trans_mag_;
	Real const rot_mag_;

	Size moved_jump_;
};

typedef utility::pointer::shared_ptr< RB_Mover > RB_MoverOP;


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

protocols::moves::TrialMoverOP
setup_MCM_trial(
	protocols::moves::MoverOP perturb,
	protocols::moves::MoverOP pack,
	protocols::moves::MoverOP minimize,
	protocols::moves::MonteCarloOP mc
);

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief mover implementation


class ProteinDNA_Relax : public protocols::moves::Mover {

public:
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::Real Real;

public:

	ProteinDNA_Relax( ScoreFunction const & scorefxn_in, Size const moving_jump_in );

	ProteinDNA_Relax( ScoreFunction const & scorefxn_in, utility::vector1< Size > const & moving_jumps );

	void
	apply_default_settings()
	{
		min_tol_ = 0.001;
		trans_mag_ = 0.5;
		rot_mag_ = 3.0;
		energycut_ = 0.1;
		ramping_initial_weight_ = 0.1;
		inner_cycles_ = 10;
		outer_cycles_ = 6;
		ramping_cycles_ = 2;
	}

	ProteinDNA_Relax &
	rot_mag( Real const setting )
	{
		rot_mag_ = setting;
		return *this;
	}

	ProteinDNA_Relax &
	trans_mag( Real const setting )
	{
		trans_mag_ = setting;
		return *this;
	}


	void
	apply( core::pose::Pose & pose );
	std::string get_name() const;

private:
	ScoreFunctionOP scorefxn_;

	utility::vector1< Size > moving_jumps_;

	Real   min_tol_; // minimization tolerance: smaller --> longer minimization
	Real trans_mag_; // controls magnitude of translation in random rigid-body perturbations
	Real   rot_mag_; // ---------------------- rotation  ---------------------------------
	Real energycut_; // controls how many positions to do rotamer trials packing

	Real ramping_initial_weight_;

	Size inner_cycles_;
	Size outer_cycles_;
	Size ramping_cycles_; // this many cycles at lower weight

};


} // namespace dna
} // namespace devel

#endif
