// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/HolesEnergy.hh
/// @brief  Packing Score
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_packing_HolesEnergy_hh
#define INCLUDED_core_scoring_packing_HolesEnergy_hh

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
// #include <core/scoring/ResidualDipolarCoupling.hh>

#include <core/scoring/packing/HolesParams.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


//Objexx headers


// Utility headers


namespace core {
namespace energy_methods {


class HolesEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:

	HolesEnergy();

	//clone
	core::scoring::methods::EnergyMethodOP
	clone() const override {
		return utility::pointer::make_shared< HolesEnergy >();
	}


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	void
	setup_for_derivatives(
		pose::Pose &,
		core::scoring::ScoreFunction const &
	)
	const override;

	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const &,
		Vector &,// F1,
		Vector & // F2
	) const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:

	core::scoring::packing::HolesParams min_params_, decoy_params_, resl_params_;
	core::Size version() const override;

};

} //scoring
} //core

#endif // INCLUDED_core_scoring_packing_HolesEnergy_HH
