// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Relax Baseclass, ClassicRelax Stage 1,2,3, ClassicRelax
///
///
/// @author Mike Tyka


#ifndef INCLUDED_protocols_relax_ClassicRelax_hh
#define INCLUDED_protocols_relax_ClassicRelax_hh


#include <protocols/relax/RelaxProtocolBase.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/RotamerTrialsMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/moves/MoverCreator.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace relax {


class ClassicRelaxCreator : public protocols::moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static std::string mover_name();
};

/// A functor class which implements the classic Rosetta++ pose_relax protocol
///
/// @todo   Crank and Wobble Moves - however there's little evidence that they are essential
///
//  @todo   Ramping of vdw --> how about providing a general ramper class. It could provide an iterator thru which
//          continuously changing score functions are accessible.
/// @brief A functor class which implements the classic Rosetta++ pose_relax protocol
class ClassicRelax: public RelaxProtocolBase {
public:
	typedef RelaxProtocolBase parent;

public:
	ClassicRelax( core::scoring::ScoreFunctionOP scorefxn_in );

	ClassicRelax( core::scoring::ScoreFunctionOP scorefxn_in, core::kinematics::MoveMapOP movemap );

	ClassicRelax( ClassicRelax const & );

	ClassicRelax();

	~ClassicRelax() override;

	protocols::moves::MoverOP clone() const override;

	static void register_options();
	// Default settings

	void set_default( core::scoring::ScoreFunctionOP scorefxn_in );

	///////////////////////////////////////////////////////////////////////////////////
	///
	/// Central Apply function
	///
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;


	///////////////////////////////////////////////////////////////////////////////////
	///
	/// Set default options from outside


	void set_default( bool const use_default_movemap = true );


	void set_default_minimizer();

	void set_default_moveset_phase1();

	void set_default_moveset_phase2();

	void set_default_moveset_phase3();

	///////////////////////////////////////////////////////////////////////////////////
	///
	/// Set options from outside


	void set_lj_ramp_cycles( int param )         { lj_ramp_cycles       = param; }

	void set_lj_ramp_inner_cycles( int param )   { lj_ramp_inner_cycles = param; }

	void set_start_rep_weight( core::Real param ){start_rep_weight      = param; }

	void set_stage2_cycles( int cycles2 )        { stage2_cycles        = cycles2; }

	void set_stage2_repack_period( int repack2 ) { stage2_repack_period = repack2; }

	void set_stage3_cycles( int cycles3 )        { stage3_cycles        = cycles3; }

	void set_tolerance( core::Real new_tolerance );

	void set_mc ( moves::MonteCarloOP new_mc_ );

	void set_full_repack ( protocols::simple_moves::PackRotamersMoverOP new_pack_full_repack );

	void set_rottrial ( protocols::simple_moves::RotamerTrialsMoverOP new_pack_rottrial );

	void setPoseExtraScore( core::pose::Pose &pose );

	void use_coarse_vdw() { //we need to know which scoring terms to ramp!
		// one could generalize this by giving a RamperClass to the Protocol, which provides an iterator that yields a scoring function
		st_rep_ = core::scoring::coarse_fa_rep;
		st_atr_ = core::scoring::coarse_fa_atr;
		st_sol_ = core::scoring::coarse_fa_sol;
	};


	///////////////////////////////////////////////////////////////////////////////////
	///
	/// Accessors


	moves::MonteCarloOP get_mc( core::pose::Pose &pose );


	protocols::checkpoint::CheckPointer & get_checkpoints() { return checkpoints_; };

private:
	// protocol stuff

	protocols::simple_moves::MinMoverOP min_mover_;
	protocols::checkpoint::CheckPointer checkpoints_;

	// THese three are special in the sense that they require a Pose at creating time.
	// this means they cannot have default values at construction time of ClassicRelax
	// because the pose is not yet known. Thus their initialisation has to be delayed
	// until apply is called. Each has a boolean flag to indicate whether the user has
	// overridden them with their own object instance or if apply should create a default
	// version.

	moves::MonteCarloOP mc_;
	bool use_default_mc_;
	void check_default_mc( core::pose::Pose & pose  );
	// stuff to do with packing

	protocols::simple_moves::PackRotamersMoverOP pack_full_repack_;
	bool use_default_pack_full_repack_;
	void check_default_full_repacker( core::pose::Pose & pose, core::kinematics::MoveMap & movemap );

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial_;
	bool use_default_pack_rottrial_;
	void check_default_rottrial( core::pose::Pose & pose, core::kinematics::MoveMap & movemap );

	// default temperature for monte carlo
	core::Real m_Temperature;

	// default number of moves
	core::Size nmoves_;
	core::Real energycut;

	// minimization options
	std::string min_type;
	bool nb_list;
	core::Real min_tolerance;

	moves::MoverOP moveset_phase1_;
	moves::MoverOP moveset_phase2_;
	moves::MoverOP moveset_phase3_;

	// PHASE1 stuff

	int lj_ramp_cycles;
	int lj_ramp_inner_cycles;
	core::Real start_rep_weight;
	core::Real end_rep_weight;

	core::scoring::ScoreType st_rep_;
	core::scoring::ScoreType st_atr_;
	core::scoring::ScoreType st_sol_;

	// PHASE2 stuff
	int stage2_repack_period;
	int stage2_cycles;

	// PHASE3 stuff RandomMover
	int stage3_cycles;


	// filters

	float score_stage2_beginning;
	float score_stage2_quarter;
	float score_stage2_half;
	float score_stage2_end;

	float filter_stage2_beginning;
	float filter_stage2_quarter;
	float filter_stage2_half;
	float filter_stage2_end;

};


}
} // protocols

#endif
