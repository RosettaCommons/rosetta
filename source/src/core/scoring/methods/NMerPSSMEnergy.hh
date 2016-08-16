// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NMerPSSMEnergy.hh
/// @brief  PSSM energy method declaration
/// @author Chris King (dr.chris.king@gmail.com)


#ifndef INCLUDED_core_scoring_methods_NMerPSSMEnergy_hh
#define INCLUDED_core_scoring_methods_NMerPSSMEnergy_hh

// Unit headers
#include <core/scoring/methods/NMerPSSMEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>
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

class NMerPSSMEnergy : public ContextIndependentOneBodyEnergy
{
public:
	typedef ContextIndependentOneBodyEnergy parent;

public:


	NMerPSSMEnergy();


	NMerPSSMEnergy( utility::vector1< std::map< core::chemical::AA, utility::vector1< core::Real > > > const & all_nmer_pssms_in );


	virtual ~NMerPSSMEnergy();

	virtual
	EnergyMethodOP
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

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }


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

	void read_nmer_pssm_list( std::string const );
	void read_nmer_pssm( std::string const );
	void nmer_length( core::Size const );
	void gate_pssm_scores( bool const );
	void nmer_pssm_scorecut( core::Real const );
	core::Real pssm_energy_at_frame_seqpos( core::Size const, core::chemical::AA const, core::Size const ) const;
	core::Size n_pssms() const;

private:
	utility::vector1< std::map< core::chemical::AA, utility::vector1< core::Real > > > all_nmer_pssms_;
	core::Size nmer_length_;
	core::Size nmer_cterm_;
	bool gate_pssm_scores_;
	core::Real nmer_pssm_scorecut_;

	void read_nmer_pssms_from_options();
	void initialize_from_options();
	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core


#endif
