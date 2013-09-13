// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/StackElecEnergy.hh
/// @brief  FA_Elec, but just 'perpendicular' to bases... trying to separate from H-bonds & geom_sol.
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_StackElecEnergy_HH
#define INCLUDED_core_scoring_rna_StackElecEnergy_HH

// Unit Headers
#include <core/scoring/rna/StackElecEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>



namespace core {
namespace scoring {
namespace rna {

typedef  numeric::xyzMatrix< Real > Matrix;

///

class StackElecEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:

	///
	StackElecEnergy( methods::EnergyMethodOptions const & options );

	///
	StackElecEnergy( StackElecEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

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


 	virtual
 	void
 	eval_atom_derivative(
 		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

  virtual
  void
  setup_for_packing(
    pose::Pose  & pose,
    utility::vector1< bool > const &,
    utility::vector1< bool > const & designing_residues
  ) const;
  
  
	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

protected:

	inline
	etable::coulomb::Coulomb const &
	coulomb() const {return coulomb_; }

private:

	//	Real
	//	get_fa_stack_score( Distance const dist, Real const cos_kappa ) const;

	Real
	get_stack_elec_score( Vector const & r_i,
												Vector const & r_j,
												Real const & i_charge,
												Real const & j_charge,
												Matrix const & M_i,
												Real & cos_kappa2 ) const;



	Vector
	get_stack_elec_deriv( Vector const & r_i,
												Vector const & r_j,
												Real const & i_charge,
												Real const & j_charge,
												Matrix const & M_i ) const;

	Real
	residue_pair_energy_one_way(
															conformation::Residue const & rsd1,
															conformation::Residue const & rsd2,
															pose::Pose const & pose,
															Real & score_base_base,
															Real & score_base_bb
															) const;

	bool
	check_base_base_OK(
	   conformation::Residue const & rsd1,
		 conformation::Residue const & rsd2,
		 Size const & m, Size const & n ) const;

	bool
	is_rna_base(
							conformation::Residue const & rsd1,
							Size const & m ) const;

	virtual
	core::Size version() const;

	etable::coulomb::Coulomb coulomb_;

	bool const base_base_only_;

	bool const verbose_;
  
  mutable bool might_be_designing_;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
