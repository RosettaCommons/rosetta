// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_etable_EtableEnergy_hh
#define INCLUDED_core_scoring_etable_EtableEnergy_hh

// Unit headers
#include <core/scoring/etable/EtableEnergy.fwd.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>


// Package headers
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/EnergyMap.hh>


namespace core {
namespace scoring {
namespace etable {


///
class EtableEnergy : public BaseEtableEnergy< EtableEnergy >  {
public:
	typedef BaseEtableEnergy< EtableEnergy >  parent;

public:
	using parent::eval_residue_pair_derivatives;
	using parent::eval_intrares_derivatives;

public:

	/// @brief construction with an etable
	EtableEnergy(
		Etable const & etable_in,
		methods::EnergyMethodOptions const & options
	);

	/// clone
	EnergyMethodOP
	clone() const {
		return new EtableEnergy( *this );
	}


	void
	setup_for_scoring_(pose::Pose const& pose, scoring::ScoreFunction const&) const;


	//
	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentTwoBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

public:

	inline
	void
	derived_prepare_for_residue_pair(
		Size const res1,
		Size const res2,
		pose::Pose const &
	) const;

private:
	// bookkeeping variable -- are intra- or interresidue weights being using at the moment?
	// don't bother to reset the BaseEtable's st_atr_, st_rep_, st_sol_ if possible.
	mutable bool using_interres_scoretypes_;
	
	virtual
	core::Size version() const;


};

/*class EtableEnergyEvaluator
{
public:
	// default ctor -- set inter-residue score types
	EtableEnergyEvaluator( EtableEnergy const & etable_energy ) :
		etable_energy_( etable_energy ),
		st_atr_( fa_atr ),
		st_rep_( fa_rep ),
		st_sol_( fa_sol )
	{}

	// switch to intra-res mode
	void set_intrares() {
		st_atr_ = fa_intra_atr;
		st_rep_ = fa_intra_rep;
		st_sol_ = fa_intra_sol;
	}

public:

	Real
	hydrogen_interaction_cutoff2() const {
		return etable_energy_.hydrogen_interaction_cutoff2();
	}

	void
	pair_energy_H(
		Atom const & at1,
		Atom const & at2,
		Real const & cp_weight,
		EnergyMap & emap
	) const;

	void
	atom_pair_energy(
		Atom const & at1,
		Atom const & at2,
		Real const & cp_weight,
		EnergyMap & emap,
		Real & dsq
	) const;


private:
	EtableEnergy const & etable_energy_;
	ScoreType st_atr_;
	ScoreType st_rep_;
	ScoreType st_sol_;
};*/

///////////////////////////////////////////////////////////////////////////////
// inline methods
///////////////////////////////////////////////////////////////////////////////

/// @brief decide between intra and interresidue score terms.  Pseudo polymorphic call
/// from base class
inline
void
EtableEnergy::derived_prepare_for_residue_pair(
	Size const res1,
	Size const res2,
	pose::Pose const &
) const {

	if ( res1 == res2 ) {
		if ( using_interres_scoretypes_ ) {
			//std::cout << "setting intaresidue scoretypes for residue " << res1 << std::endl;
			set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
			using_interres_scoretypes_ = false;
		}
		assert( rep_scoretype() == fa_intra_rep );
	} else {
		if ( ! using_interres_scoretypes_ ) {
			//std::cout << "setting interresidue scoretypes for residue pair " << res1 << " & " << res2 << std::endl;
			set_scoretypes( fa_atr, fa_rep, fa_sol );
			using_interres_scoretypes_ = true;
		}
		assert( rep_scoretype() == fa_rep );
	}

}

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
