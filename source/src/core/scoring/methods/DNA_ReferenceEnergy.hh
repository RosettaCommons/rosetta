// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/methods/DNA_ReferenceEnergy.hh
/// @brief  dna scoring
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_DNA_ReferenceEnergy_HH
#define INCLUDED_core_scoring_methods_DNA_ReferenceEnergy_HH

// Unit Headers
#include <core/scoring/methods/DNA_ReferenceEnergy.fwd.hh>
#include <core/scoring/dna/BasePartner.hh> // inline fxn


// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
// #include <core/pack/task/PackerTask.hh>
//#include <core/scoring/dna/DNA_ReferencePotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>

// Utility headers


namespace core {
namespace scoring {
namespace methods {

///
class DNA_ReferenceEnergy : public ContextIndependentTwoBodyEnergy {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;
public:

	///
	DNA_ReferenceEnergy( EnergyMethodOptions const & options );

	///
	DNA_ReferenceEnergy( DNA_ReferenceEnergy const & src );


	/// clone
	virtual
	EnergyMethodOP
	clone() const;


	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	// virtual
	// void
	// setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	// virtual
	// void
	// setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	// virtual
	// void
	// setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////


	/// helper functions, possible used outside also
	inline
	Real
	base_step_energy( chemical::AA const & first_aa, chemical::AA const & second_aa ) const
	{
		return base_step_reference_energies_
			[ int(  first_aa ) - chemical::first_DNA_aa + 1 ]
			[ int( second_aa ) - chemical::first_DNA_aa + 1 ];
	}

	/// ordered!!!! requires pos1<pos2
	inline
	bool
	count_pair_bs(
		Size const pos1,
		Size const pos2,
		core::scoring::dna::BasePartner const & partner
	) const
	{
		return ( pos2 == pos1 + 1 &&
			pos1 >= 1 && pos2 <= partner.size() &&
			partner[pos1] && partner[pos2] &&
			partner[pos2] == partner[pos1]-1 &&
			partner[pos1] != pos2 );
	}


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
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}


	// virtual
	// void
	// eval_atom_derivative(
	//  id::AtomID const & atom_id,
	//  pose::Pose const & pose,
	//  kinematics::DomainMap const & domain_map,
	//  ScoreFunction const &,
	//  EnergyMap const & weights,
	//  Vector & F1,
	//  Vector & F2
	// ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	///@brief Return the version of the energy method
	virtual
	core::Size version() const { return 1; }


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	utility::vector1< utility::vector1< Real > > base_step_reference_energies_;

	// bool use_water_anchor_atom_weights_;
	// std::map< std::string, Real > water_anchor_atom_weights_;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
