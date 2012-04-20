// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/NMerRefEnergy.hh
/// @brief  Reference energy method declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_NMerRefEnergy_hh
#define INCLUDED_core_scoring_methods_NMerRefEnergy_hh

// Unit headers
#include <core/scoring/methods/NMerRefEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <map>


// Utility headers


namespace core {
namespace scoring {
namespace methods {

class NMerRefEnergy : public ContextIndependentOneBodyEnergy
{
public:
	typedef ContextIndependentOneBodyEnergy parent;

public:

	///
	NMerRefEnergy();

	///
	NMerRefEnergy( std::map< std::string, core::Real > const & nmer_ref_energies_in );

	///
	virtual ~NMerRefEnergy();

	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	///
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	///
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


private:
	std::map< std::string, core::Real > nmer_ref_energies_;
	core::Size nmer_length_;
	core::Size nmer_nterm_;
	core::Size nmer_cterm_;

	void read_nmer_energy_table();
	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core


#endif
