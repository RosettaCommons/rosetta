// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CSD_TorsionEnergy.hh
///
/// @brief
/// @author Ian W. Davis
/// @author Kristian Kaufmann


#ifndef INCLUDED_core_scoring_methods_CSD_TorsionEnergy_hh
#define INCLUDED_core_scoring_methods_CSD_TorsionEnergy_hh

#include <core/scoring/methods/CSD_TorsionEnergy.fwd.hh>

#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

namespace core {
namespace scoring {
namespace methods {


///@brief A knowledge-based torsional potential for small molecules
/// derived from the Cambridge Structural Database by KWK.
///
///@details Implemented as a two-body term to allow for torsions
/// between multi-residue ligands.
///
class CSD_TorsionEnergy : public ContextIndependentTwoBodyEnergy
{
public:

	CSD_TorsionEnergy();
	virtual ~CSD_TorsionEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/// setup for packing
	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	/// setup for scoring
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	/// setup for derivatives
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & ) const ;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// @brief CSD_TorsionEnergy does not have an atomic interation threshold
	virtual
	Distance
	atomic_interaction_cutoff() const;

	/// @brief CSD_TorsionEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

private:


}; // CSD_TorsionEnergy


} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_CSD_TorsionEnergy_HH
