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
/// @author Liz Kellogg

#ifndef INCLUDED_protocols_simple_moves_symmetry_SymShakeStructureMover_hh
#define INCLUDED_protocols_simple_moves_symmetry_SymShakeStructureMover_hh

#include <protocols/symmetry/SymShakeStructureMover.fwd.hh>
#include <protocols/simple_moves/ShakeStructureMover.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

// C++ headers

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

//protocols
#include <protocols/moves/Mover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace symmetry {

class SymShakeStructureMover: public protocols::simple_moves::ShakeStructureMover{
public:
	SymShakeStructureMover();

	SymShakeStructureMover(core::scoring::ScoreFunctionOP s);

	SymShakeStructureMover(core::scoring::ScoreFunctionOP s,
		core::Real temperature);

	SymShakeStructureMover(core::scoring::ScoreFunctionOP s,
		core::Real ens_diversity,
		core::Real ens_div_tolerance);

	virtual ~SymShakeStructureMover();

	virtual std::string get_name() const;
	//setters
	//   void set_scorefunction(core::scoring::ScoreFunction & s);

	void
	reduce_fa_rep(float fraction_fa_rep, core::scoring::ScoreFunction & s);

	void minimize_with_constraints(core::pose::Pose & p,
		core::scoring::ScoreFunction & s);

	void setup_for_run(core::pose::Pose & p);

	void run_mc(core::pose::Pose & p, core::scoring::ScoreFunction & s,
		core::Real temperature);

};

} //symmetry
} //protocols

#endif
