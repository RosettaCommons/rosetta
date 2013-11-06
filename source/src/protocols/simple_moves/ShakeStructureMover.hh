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
/// @author Liz Kellogg

#ifndef INCLUDED_protocols_simple_moves_ShakeStructureMover_HH
#define INCLUDED_protocols_simple_moves_ShakeStructureMover_HH

#include <protocols/simple_moves/ShakeStructureMover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

// C++ headers
// AUTO-REMOVED #include <cstdlib>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//protocols
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

	// no using in .hh files !
	// using namespace core;
	// using namespace scoring;

class ShakeStructureMover: public protocols::moves::Mover{
public:
	ShakeStructureMover();

	ShakeStructureMover(core::scoring::ScoreFunctionOP s);

	ShakeStructureMover(core::scoring::ScoreFunctionOP s,
		core::Real temperature);

	ShakeStructureMover(core::scoring::ScoreFunctionOP s,
		core::Real ens_diversity,
		core::Real ens_div_tolerance);

	virtual ~ShakeStructureMover();

	//setters
	void set_skip_low_temp_phase( bool truefalse);
	void set_mc_temperature(core::Real temp);
	void set_nrounds( int new_nrounds );
	void set_ramp_fa_rep(bool truefalse);
	void set_minimize_with_cst(bool truefalse);
	void set_scorefunction(core::scoring::ScoreFunction & s);
	void set_ensemble_diversity(core::Real ca_rmsd);
	void set_rmsd_target_tolerance(core::Real tol);
	void set_sc_min(bool truefalse);

	//accessors
	core::Real get_mc_temperature();
	bool get_ramp_fa_rep();
	bool get_minimize_with_cst();
	core::scoring::ScoreFunctionOP get_scorefunction();
	core::Real get_ensemble_diversity();
	core::Real get_rmsd_target_tolerance();
	bool get_sc_min();
	//
	core::Real get_harmonic_ca_cst_std_dev();
	core::Real get_ensemble_ca_rmsd();
	bool get_skip_low_temp_phase();
	bool get_min_cst();
	bool get_testing_phase();
	bool get_scorefunction_initialized();
	core::scoring::ScoreFunctionOP get_min_scorefunction();
	core::Size get_nrounds();
	//run-time

	void apply(core::pose::Pose & p);
	virtual std::string get_name() const;

protected:
	// for derived classes
	void set_testing_phase( bool truefalse );
	void set_mc_temp( core::Real temperature );
	void set_is_properly_initialized( bool truefalse );
	void set_min_scorefunction( core::scoring::ScoreFunctionOP scfxn );

private:
	core::Real mc_temp;
	bool ramp_fa_rep;
	bool min_cst;
	core::scoring::ScoreFunctionOP scorefxn;
	core::Real ensemble_ca_rmsd;
	core::Real ensemble_ca_rmsd_tolerance;
	bool is_properly_initialized;
	core::Real harmonic_ca_cst_std_dev;
	bool scorefunction_initialized;
	bool sc_min_only;
	int nrounds;
	double cst_weight;
	bool skip_low_temp_phase;
	core::scoring::ScoreFunctionOP min_scorefxn; //to avoid having to reinitialize
	//every time we minimize
	bool testing_phase;

protected:
		void
		reduce_fa_rep(float fraction_fa_rep, core::scoring::ScoreFunction & s);

	virtual void
		setup_for_run(core::pose::Pose & p);

	virtual void
		minimize_with_constraints(core::pose::Pose & p,
		core::scoring::ScoreFunction & s);

	virtual void
		setup_ca_constraints(core::pose::Pose & pose,
		core::scoring::ScoreFunction & s,
		float const CA_cutoff,
		float const cst_tol);

	virtual void
		run_mc(core::pose::Pose & p, core::scoring::ScoreFunction & s,
		core::Real temperature);

	core::Real
		set_temp_based_on_ens_diversity(core::pose::Pose & p,
		core::scoring::ScoreFunction & s);

	void
	setup_movers(
		protocols::simple_moves::SmallMoverOP gsmall,
		protocols::simple_moves::ShearMoverOP gshear,
		core::Real small_H_angle_max, core::Real small_E_angle_max, core::Real small_L_angle_max,
		core::Real shear_H_angle_max, core::Real shear_E_angle_max, core::Real shear_L_angle_max);

	void
		create_ensemble(core::pose::Pose & p,
		core::scoring::ScoreFunction & s);

};

} //protocols
} //simple_moves

#endif //INCLUDED_protocols_simple_moves_ShakeStructureMover_HH
