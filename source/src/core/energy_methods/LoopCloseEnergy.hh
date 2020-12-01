// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LoopCloseEnergy.hh
/// @brief  Loop closure energy, currently defined for RNA.
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_loop_graph_LoopCloseEnergy_hh
#define INCLUDED_core_scoring_loop_graph_LoopCloseEnergy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/loop_graph/LoopGraph.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


class LoopCloseEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:

	/// @brief Defines a loop closure energy based on FullModelInfo and pose geometry.
	LoopCloseEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	/// copy ctor
	LoopCloseEnergy( LoopCloseEnergy const & src );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	/////////////////////////////////
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2 ) const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:

	void
	update_loop_atoms_and_lengths( pose::Pose & pose ) const;

	Real
	get_loop_distance2( Vector const & v_takeoff,  Vector const & v_landing ) const;

	Vector
	get_loop_distance2_deriv( Vector const & v_takeoff,  Vector const & v_landing ) const;

	core::Size version() const override;

	core::scoring::methods::EnergyMethodOptions const & options_;

	mutable core::scoring::loop_graph::LoopGraphOP loop_graph_;

};


}
}

#endif // INCLUDED_core_energy_methods_LoopCloseEnergy_HH
