// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNAChiEnergy.hh
/// @brief  Energy term for simple DNA chi torsion
/// @author Jim Havranek


#ifndef INCLUDED_core_scoring_dna_DNAChiEnergy_HH
#define INCLUDED_core_scoring_dna_DNAChiEnergy_HH

// Unit headers
#include <core/scoring/dna/DNAChiEnergy.fwd.hh>
#include <core/scoring/dna/DNABFormPotential.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>


namespace core {
namespace scoring {
namespace dna {


class DNAChiEnergy : public methods::ContextIndependentOneBodyEnergy {
public:
	typedef ContextIndependentOneBodyEnergy parent;
public:

	/// @brief ctor
	DNAChiEnergy();

	/// @brief dtor
	virtual ~DNAChiEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
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

	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;


	virtual
	core::Size version() const;

	// data
private:
	// Will probably need to store data here, although I may hard-wire this for the first pass
	core::scoring::dna::DNABFormPotential const & potential_;

};

} // dna
} // scoring
} // core


#endif // INCLUDED_core_scoring_dna_DNAChiEnergy_HH
