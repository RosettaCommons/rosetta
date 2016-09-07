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
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymFoldAndDockRbTrialMover_hh
#define INCLUDED_protocols_symmetric_docking_SymFoldAndDockRbTrialMover_hh

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockRbTrialMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>


// Utility Headers

namespace protocols {
namespace symmetric_docking {

class SymFoldandDockRbTrialMover : public moves::Mover
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

	// init to defaults
	void init();

	~SymFoldandDockRbTrialMover() override= default;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool smooth_move_;
	core::Real rot_mag_, trans_mag_;
	core::Real rot_mag_smooth_, trans_mag_smooth_;
	core::Size rigid_body_cycles_;
	bool mc_filter_, rotate_anchor_to_x_;
};

} // symmetric_docking
} // rosetta
#endif
