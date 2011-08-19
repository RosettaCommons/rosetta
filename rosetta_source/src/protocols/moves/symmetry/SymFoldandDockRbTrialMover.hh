// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_moves_symmetry_SymFoldandDockRbTrialMover_hh

#define INCLUDED_protocols_moves_symmetry_SymFoldandDockRbTrialMover_hh

// Unit headers
#include <protocols/moves/symmetry/SymFoldandDockRbTrialMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Utility Headers

namespace protocols {
namespace moves {
namespace symmetry {
///////////////////////////////////////////////////////////////////////////////
class SymFoldandDockRbTrialMover : public Mover
{
public:

typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;

public:

	// default constructor
	SymFoldandDockRbTrialMover();

	SymFoldandDockRbTrialMover( core::scoring::ScoreFunctionCOP scorefxn );

	SymFoldandDockRbTrialMover( core::scoring::ScoreFunctionCOP scorefxn, bool smooth_move );

	SymFoldandDockRbTrialMover(
		core::scoring::ScoreFunctionCOP scorefxn,
		bool smooth_move,
		core::Real rot_mag,
		core::Real trans_mag
	);

	~SymFoldandDockRbTrialMover(){}

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool smooth_move_;
	core::Real rot_mag_;
	core::Real trans_mag_;
	core::Size rigid_body_cycles_;
	bool mc_filter_;
};

} // symmetry
} // moves
} // rosetta
#endif
