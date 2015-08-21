// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/HolesEnergyRes.hh
/// @brief  Packing Score
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_packing_HolesEnergyRes_hh
#define INCLUDED_core_scoring_packing_HolesEnergyRes_hh

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
//#include <core/scoring/ResidualDipolarCoupling.hh>

#include <core/scoring/packing/HolesParamsRes.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreType.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

//Objexx headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>


// Utility headers


namespace core {
namespace scoring {
namespace packing {


class HolesEnergyRes : public methods::WholeStructureEnergy {
public:
	typedef methods::WholeStructureEnergy parent;

public:

	HolesEnergyRes();

	//clone
	virtual
	methods::EnergyMethodOP
	clone() const {
		return methods::EnergyMethodOP( new HolesEnergyRes() );
	}


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	void
	setup_for_derivatives(
		pose::Pose &,
		ScoreFunction const &
	)
	const;

	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const &,
		Vector &,// F1,
		Vector & // F2
	) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:

	HolesParamsRes params_;
	virtual
	core::Size version() const;
};

} //packing
} //scoring
} //core

#endif // INCLUDED_core_scoring_packing_HolesEnergyRes_HH
