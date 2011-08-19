// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PoissonBoltzmannEnergy.hh
/// @brief  Membrane Environment Cbeta Energy
/// @author Bjorn Wallner


#ifndef INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_HH
#define INCLUDED_core_scoring_methods_PoissonBoltzmannEnergy_HH

// Unit Headers
#include <core/scoring/methods/PoissonBoltzmannEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


namespace core {
namespace scoring {
namespace methods {

///
class PoissonBoltzmannEnergy : public ContextIndependentLRTwoBodyEnergy {
public:
	typedef ContextIndependentLRTwoBodyEnergy  parent;

public:

	///
	PoissonBoltzmannEnergy();


	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	methods::LongRangeEnergyType
	long_range_type() const;

	virtual bool defines_intrares_energy( EnergyMap const &  ) const { return true; }

	virtual bool defines_residue_pair_energy(
											 pose::Pose const & pose,
											 Size res1,
											 Size res2
											 ) const;
	
	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;


	virtual
	void
	eval_intrares_energy(
						 conformation::Residue const & rsd,
						 pose::Pose const & pose,
						 ScoreFunction const & sfxn,
						 EnergyMap & emap
						 ) const ;

	///
	virtual void residue_pair_energy(
									 conformation::Residue const & rsd1,
									 conformation::Residue const & rsd2,
									 pose::Pose const & pose,
									 ScoreFunction const & sfxn,
									 EnergyMap & emap
									 ) const;

	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	//virtual
	//void
	//finalize_total_energy(
	//	pose::Pose & pose,
	//	ScoreFunction const &,
	//	EnergyMap &
	//) const;

	/// @brief PB Energy is context independent and thus indicates that no context graphs need to
	/// be maintained by class Energies
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const;
	

/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////

private:

// const-ref to scoring database
//core::scoring::PoissonBoltzmannPotential const & potential_;
bool potential_is_loaded_;
core::Size fixed_residue_;
	
virtual
core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
