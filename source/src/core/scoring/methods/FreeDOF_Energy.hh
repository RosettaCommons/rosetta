// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FreeDOF_Energy.hh
/// @brief  Score function class
/// @author Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_scoring_methods_FreeDOF_Energy_hh
#define INCLUDED_core_scoring_methods_FreeDOF_Energy_hh

// Unit headers
#include <core/scoring/methods/FreeDOF_Energy.fwd.hh>
#include <core/scoring/methods/FreeDOF_Options.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class FreeDOF_Energy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief ctor
	FreeDOF_Energy( EnergyMethodOptions const & energy_method_options );

	/// @brief dtor
	virtual ~FreeDOF_Energy();

	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;


	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap & emap
	) const;


	/// @brief FreeDOF_Energy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	virtual
	core::Size version() const;

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals ) const;

private:


	void
	accumulate_stack_energy(
		pose::Pose & pose,
		ScoreFunction const & scorefxn,
		utility::vector1< Real > & stack_energy ) const;

	void
	do_fa_stack_scorefunction_checks( ScoreFunction const & scorefxn ) const;

	void
	get_hbond_energy(
		pose::Pose & pose,
		ScoreFunction const & scorefxn,
		utility::vector1< Real > & base_hbond_energy,
		utility::vector1< Real > & sugar_hbond_energy ) const;

private:

	EnergyMethodOptions const & energy_method_options_;
	FreeDOF_Options const & options_;
	utility::vector1< Real > const & free_res_weights_;

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
