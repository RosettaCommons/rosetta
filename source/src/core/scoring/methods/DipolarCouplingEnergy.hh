// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/DipolarCouplingEnergy.hh
/// @brief  CSP energy
/// @author Lei Shi


#ifndef INCLUDED_core_scoring_methods_DipolarCouplingEnergy_hh
#define INCLUDED_core_scoring_methods_DipolarCouplingEnergy_hh

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/DipolarCoupling.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
//#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID_Map.hh>


//Objexx headers


// Utility headers


namespace core {
namespace scoring {
namespace methods {


class DipolarCouplingEnergy : public WholeStructureEnergy  {
public:
	typedef WholeStructureEnergy  parent;

public:

	DipolarCouplingEnergy();

	//clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////
	virtual
	void
	setup_for_scoring( pose::Pose &, ScoreFunction const & ) const;

	/// @brief Called at the beginning of atom tree minimization, this method
	/// allows the derived class the opportunity to initialize pertinent data
	/// that will be used during minimization.  During minimzation, the chemical
	/// structure of the pose is constant, so assumptions on the number of atoms
	/// per residue and their identities are safe so long as the pose's Energies
	/// object's "use_nblist()" method returns true.
	virtual
	void
	setup_for_minimizing(
		pose::Pose & ,
		ScoreFunction const & ,
		kinematics::MinimizerMapBase const &
	) const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}


	virtual void eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

private:

	DipolarCoupling& dc_from_pose(
		pose::Pose & pose
	) const;

	Real eval_dc(
		pose::Pose & pose
	) const;


private:

	//used by Energy Method during scoring... should this become part of DipolarCoupling and thus cached in the pose
	mutable core::Real dc_score_; //computed in setup_for_scoring.. delivered in finalize
	mutable id::AtomID_Map< utility::vector1<Size> > atom2dc_map_;
	virtual
	core::Size version() const;
};

} //methods
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
