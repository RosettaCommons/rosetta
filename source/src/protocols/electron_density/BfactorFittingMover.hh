// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file

#ifndef INCLUDED_protocols_electron_density_BfactorFittingMover_hh
#define INCLUDED_protocols_electron_density_BfactorFittingMover_hh

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>

namespace protocols {
namespace electron_density {


/// Bfactor multifunc
class BfactorMultifunc : public core::optimization::Multifunc {
public:
	BfactorMultifunc(
		core::pose::Pose & pose_in,
		core::Real wt_adp,
		core::Real wt_dens,
		core::Real rmax,
		core::Real radius_exp,
		core::Real scorescale,
		bool exact );

	virtual ~BfactorMultifunc() {}

	virtual	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	virtual	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const;

	// helper functions
	void
	poseBfacts2multivec( core::pose::Pose & pose, core::optimization::Multivec &y ) const;

	void
	multivec2poseBfacts( core::optimization::Multivec const &vars, core::pose::Pose & pose ) const;

private:
	core::pose::Pose & pose_;

	core::Real scorescale_;  // overall weighing on score
	core::Real wt_adp_;      // weight on constraints
	core::Real wt_dens_;     // weight on density
	core::Real rmax_;        // max radius (def 5)
	core::Real radius_exp_;  // weigh neighbors by 1/r^thisval (def 1)
	bool exact_;

	// handles AtomID <-> vector index mapping
	core::id::AtomID_Map< core::Size > atom_indices_;
	utility::vector1< core::id::AtomID > moving_atoms_;

	//
	core::Real B_EPS;
};


/// mover to fit B factors
class BfactorFittingMover : public moves::Mover {
public:
	BfactorFittingMover();
	virtual ~BfactorFittingMover() {}

	void init();

	virtual void apply( core::pose::Pose & );

	std::string get_name() const { return "BfactorFitting"; }

	moves::MoverOP clone() const { return new BfactorFittingMover( *this ); }
	moves::MoverOP fresh_instance() const { return new BfactorFittingMover; }

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & );

private:
	core::Real scorescale_;  // overall weighing on score
	core::Real wt_adp_;      // weight on constraints
	core::Real wt_dens_;     // weight on density
	core::Real rmax_;        // max radius (def 5)
	core::Size radius_exp_;  // weigh neighbors by 1/r^thisval (def 1)
	core::Size max_iter_;    // minimizer iterations
	std::string minimizer_;
	bool init_;
	bool exact_;

	bool opt_to_fsc_;
	core::Real res_low_,res_high_;
	core::Size nresbins_;

	utility::vector1< core::Real > weights_to_scan_;
};

} // moves
} // protocols

#endif
