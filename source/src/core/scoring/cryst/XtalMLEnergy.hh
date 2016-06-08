// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/cryst/XtalMLEnergy.hh
/// @brief  ML target
/// @author Frank DiMaio


#ifndef INCLUDED_core_scoring_cryst_XtalMLEnergy_hh
#define INCLUDED_core_scoring_cryst_XtalMLEnergy_hh

// Package headers
#include <core/scoring/cryst/XtalMLEnergyCreator.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <string>

namespace core {
namespace scoring {
namespace cryst {

class XtalMLEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:


	XtalMLEnergy();


	XtalMLEnergy( XtalMLEnergy const & src ) :
		parent(src),
		dml_dx(src.dml_dx),
		ml(src.ml)
	{}

	virtual methods::EnergyMethodOP clone() const;

	virtual void finalize_total_energy(pose::Pose & pose,ScoreFunction const &,EnergyMap & totals) const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & scorefxn ) const;

	virtual void setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	virtual void
	setup_for_minimizing( pose::Pose & pose, ScoreFunction const & sf, kinematics::MinimizerMapBase const & ) const;

	virtual void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	virtual void indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:
	//////
	// PRIVATE DATA
	//////
	// precomputed derivatives
	mutable utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > dml_dx;

	// ml target function
	mutable core::Real ml;

	virtual
	core::Size version() const;
};

}
}
}

#endif
